import streamlit as st
import math
import pandas as pd
import io
import sqlite3
import json

# ==============================================================================
# DATABASE LOGIC (New Section)
# ==============================================================================
DB_NAME = "solar_projects.db"

def init_db():
    conn = sqlite3.connect(DB_NAME)
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS projects 
                 (name TEXT PRIMARY KEY, timestamp DATETIME DEFAULT CURRENT_TIMESTAMP, data_json TEXT)''')
    conn.commit()
    conn.close()

def save_project_to_db(name, params):
    # Convert NaN to None for JSON compatibility
    clean_params = params.copy()
    for k, v in clean_params.items():
        if isinstance(v, float) and math.isnan(v):
            clean_params[k] = None
            
    conn = sqlite3.connect(DB_NAME)
    c = conn.cursor()
    data_str = json.dumps(clean_params)
    c.execute("INSERT OR REPLACE INTO projects (name, data_json) VALUES (?, ?)", (name, data_str))
    conn.commit()
    conn.close()

def load_project_from_db(name):
    conn = sqlite3.connect(DB_NAME)
    c = conn.cursor()
    c.execute("SELECT data_json FROM projects WHERE name = ?", (name,))
    row = c.fetchone()
    conn.close()
    if row:
        data = json.loads(row[0])
        # Convert None back to NaN
        for k, v in data.items():
            if v is None:
                data[k] = float('nan')
        return data
    return None

def list_projects():
    conn = sqlite3.connect(DB_NAME)
    c = conn.cursor()
    c.execute("SELECT name FROM projects ORDER BY timestamp DESC")
    names = [row[0] for row in c.fetchall()]
    conn.close()
    return names

# Initialize DB
init_db()

# ==============================================================================
# ORIGINAL LOGIC (No changes made to logic below)
# ==============================================================================

# Define Euler's Constant
EULER_E = 2.718281828

# --- Define Lookup Table for a, b, DeltaTcnd ---
COEFFICIENT_LOOKUP = {
    ("Glass/cell/glass", "Open rack"): {"a": -3.47, "b": -0.0594, "delta_tcnd": 3},
    ("Glass/cell/glass", "Close-roof mount"): {"a": -2.98, "b": -0.0471, "delta_tcnd": 1},
    ("Glass/cell/polymer sheet", "Open rack"): {"a": -3.56, "b": -0.0750, "delta_tcnd": 3},
    ("Glass/cell/polymer sheet", "Insulated back"): {"a": -2.81, "b": -0.0455, "delta_tcnd": 0},
    ("Polymer/thin-film/steel", "Open rack"): {"a": -3.58, "b": -0.1130, "delta_tcnd": 3},
}

# --- Helper function for conditional formatting of values ---
def format_value(value, decimal_places=0):
    if math.isinf(value):
        return "Inf"
    if math.isnan(value):
        return "N/A"
    return f"{value:.{decimal_places}f}"

# --- Helper function to parse comma-separated numbers ---
def parse_numbers_from_string(input_string, default_value):
    if not input_string:
        return [default_value]
    numbers = []
    parts = input_string.split(',')
    for part in parts:
        try:
            num = int(part.strip())
            numbers.append(num)
        except ValueError:
            st.warning(f"Skipping invalid input: '{part.strip()}' is not a valid number. Using default {default_value}.")
    return numbers if numbers else [default_value]

# --- Helper function to calculate Tm and Tcell for a single row ---
def calculate_tcell_for_row_internal(row, a, b, delta_tcnd, gain):
    ws_val = row['WS']
    ws_adj = max(ws_val, 0.1) if ws_val is not None else 0.1

    ghi_val = row['GHI'] if row['GHI'] is not None else 0

    exp_term_arg = a + b * ws_adj
    if exp_term_arg > 700:
        exp_term = float('inf')
    elif exp_term_arg < -700:
        exp_term = 0.0
    else:
        exp_term = EULER_E**(exp_term_arg)

    tm_val = (ghi_val * gain * exp_term) + row['TEMP']
    tcell_val = tm_val + ((ghi_val * gain) / 1000) * delta_tcnd

    return tcell_val

def calculate_solar_pv_design(
    module_supplier, module_type, module_vmpp, module_voc, module_impp, module_isc, module_power_stc, module_v_max_system,
    module_temp_coeff_pmax, module_temp_coeff_voc, module_temp_coeff_isc, module_noct,
    module_dim_width, module_dim_length,
    selected_coeff_a, selected_coeff_b, selected_coeff_delta_tcnd,
    inverter_supplier, inverter_type, inverter_transformer_integrated,
    inverter_vmpp_min, inverter_vmpp_max, inverter_v_system_max,
    inverter_max_recommended_pv_power_kw, inverter_nominal_pv_power_kw,
    inverter_max_pv_current_a, inverter_nominal_pv_current_a, inverter_nb_inputs_cc, inverter_isc_max_per_inputs,
    design_azimuth, design_tilt_angle, design_row_spacing_m,
    design_pv_module_rated_power_wp,
    design_modules_per_string, design_strings_per_inverter, design_num_inverters, design_inverter_rated_ac_power_kVA,
    max_op_temp_c, min_op_temp_c,
    max_temp_inclination_gain, min_temp_inclination_gain
):
    results = {}
    STC_temp = 25

    results['Design Modules per String'] = design_modules_per_string
    results['Design Strings per Inverter'] = design_strings_per_inverter

    results['Min temp for Voc to reach max inverter voltage (°C)'] = float('nan')
    results['Min temp for Vmpp to reach MPPT limit (upper) (°C)'] = float('nan')
    results['Max temp for Vmpp to reach MPPT limit (lower) (°C)'] = float('nan')
    results['Max temp for Isc to reach limit (°C)'] = float('nan')

    total_rated_power_PDC_Wp = (
        design_modules_per_string * design_strings_per_inverter * design_num_inverters * design_pv_module_rated_power_wp
    )
    results['Total rated power PDC (MWp)'] = total_rated_power_PDC_Wp / 1_000_000

    total_rated_power_PAC_MW = design_inverter_rated_ac_power_kVA / 1000
    results['Total rated power PAC (MW)'] = total_rated_power_PAC_MW

    if total_rated_power_PAC_MW != 0:
        results['Ratio PDC PAC'] = results['Total rated power PDC (MWp)'] / total_rated_power_PAC_MW
    else:
        results['Ratio PDC PAC'] = 0

    theoretical_max_modules_voc = inverter_v_system_max / module_voc
    rounded_max_modules_voc = round(theoretical_max_modules_voc)
    if rounded_max_modules_voc * module_voc > inverter_v_system_max:
        max_modules_voc_calc = rounded_max_modules_voc - 1
    else:
        max_modules_voc_calc = rounded_max_modules_voc
    results['Max string (Voc_limit_calc)'] = max_modules_voc_calc

    theoretical_min_modules_vmpp = inverter_vmpp_min / module_vmpp
    rounded_min_modules_vmpp = round(theoretical_min_modules_vmpp)
    if rounded_min_modules_vmpp * module_vmpp < inverter_vmpp_min:
        min_modules_vmpp_calc = rounded_min_modules_vmpp + 1
    else:
        min_modules_vmpp_calc = rounded_min_modules_vmpp
    results['Min string (Vmpp_min_limit_calc)'] = min_modules_vmpp_calc

    theoretical_max_modules_vmpp = inverter_vmpp_max / module_vmpp
    rounded_max_modules_vmpp = round(theoretical_max_modules_vmpp)
    if rounded_max_modules_vmpp * module_vmpp > inverter_vmpp_max:
        max_modules_vmpp_calc = rounded_max_modules_vmpp - 1
    else:
        max_modules_vmpp_calc = rounded_max_modules_vmpp
    results['Max string (Vmpp_max_limit_calc)'] = max_modules_vmpp_calc

    if module_impp != 0:
        max_strings_per_mppt_limit_calc = math.floor(inverter_max_pv_current_a / module_impp)
    else:
        max_strings_per_mppt_limit_calc = 0
    results['Max string (p/MPPT)_limit_calc'] = max_strings_per_mppt_limit_calc

    results['Calculated String Voc (STC) for Max String (N16)'] = max_modules_voc_calc * module_voc
    results['Check Voc (STC) vs Inverter Max (N17)'] = "OK" if results['Calculated String Voc (STC) for Max String (N16)'] < inverter_v_system_max else "NOK"

    results['Calculated String Vmpp (STC) for Min String (O16)'] = min_modules_vmpp_calc * module_vmpp
    results['Check Vmpp (STC) vs Inverter Min (O17)'] = "OK" if results['Calculated String Vmpp (STC) for Min String (O16)'] > inverter_vmpp_min else "NOK"

    results['Calculated String Vmpp (STC) for Max String (P16)'] = max_modules_vmpp_calc * module_vmpp
    results['Check Vmpp (STC) vs Inverter Max (P17)'] = "OK" if results['Calculated String Vmpp (STC) for Max String (P16)'] < inverter_vmpp_max else "NOK"

    results['Configured_num_strings_per_inverter'] = design_strings_per_inverter
    dc_power_per_inverter_kWp = (design_modules_per_string * design_strings_per_inverter * module_power_stc) / 1000
    results['Configured_DC_Power_per_inverter_kWp'] = dc_power_per_inverter_kWp

    results['Check_Config_DC_Power'] = "OK" if dc_power_per_inverter_kWp <= inverter_max_recommended_pv_power_kw else "NOK"

    temp_diff_max_op = max_op_temp_c - STC_temp
    voc_temp_factor_max = 1 + (module_temp_coeff_voc * temp_diff_max_op / 100)
    string_voc_at_max_op_temp = design_modules_per_string * module_voc * voc_temp_factor_max
    results['String Voc (V) at Max Op Temp'] = string_voc_at_max_op_temp

    vmpp_temp_factor_max = 1 + (module_temp_coeff_voc * temp_diff_max_op / 100)
    string_vmpp_at_max_op_temp = design_modules_per_string * module_vmpp * vmpp_temp_factor_max
    results['String Vmpp (V) at Max Op Temp'] = string_vmpp_at_max_op_temp

    isc_temp_factor_max = 1 + (module_temp_coeff_isc * temp_diff_max_op / 100)
    array_isc_at_max_op_temp = (
        design_strings_per_inverter * module_isc * isc_temp_factor_max
    )
    results['Array Isc (A) at Max Op Temp'] = array_isc_at_max_op_temp

    temp_diff_min_op = min_op_temp_c - STC_temp
    voc_temp_factor_min = 1 + (module_temp_coeff_voc * temp_diff_min_op / 100)
    string_voc_at_min_op_temp = design_modules_per_string * module_voc * voc_temp_factor_min
    results['String Voc (V) at Min Op Temp'] = string_voc_at_min_op_temp

    vmpp_temp_factor_min = 1 + (module_temp_coeff_voc * temp_diff_min_op / 100)
    string_vmpp_at_min_op_temp = design_modules_per_string * module_vmpp * vmpp_temp_factor_min
    results['String Vmpp (V) at Min Op Temp'] = string_vmpp_at_min_op_temp

    isc_temp_factor_min = 1 + (module_temp_coeff_isc * temp_diff_min_op / 100)
    array_isc_at_min_op_temp = (
        design_strings_per_inverter * module_isc * isc_temp_factor_min
    )
    results['Array Isc (A) at Min Op Temp'] = array_isc_at_min_op_temp

    results['string_voc_change_per_degC'] = (module_temp_coeff_voc / 100) * module_voc * design_modules_per_string
    results['temp_diff_max_op_exact'] = max_op_temp_c - STC_temp
    results['total_voc_change_max_temp'] = results['temp_diff_max_op_exact'] * results['string_voc_change_per_degC']
    results['string_voc_stc'] = module_voc * design_modules_per_string
    results['array_current_change_per_degC'] = (module_temp_coeff_isc / 100) * module_impp * design_strings_per_inverter
    results['array_impp_stc'] = module_impp * design_strings_per_inverter
    results['total_current_change_max_temp'] = inverter_max_pv_current_a - results['array_impp_stc']
    results['array_impp_at_max_op_temp'] = results['array_impp_stc'] + (results['temp_diff_max_op_exact'] * results['array_current_change_per_degC'])

    temp_limit_voc_max_inv_val = float('nan')
    if module_temp_coeff_voc != 0:
        temp_limit_voc_max_inv_val = (
            (100 * ((inverter_v_system_max / (design_modules_per_string * module_voc)) - 1))
            / module_temp_coeff_voc
        ) + STC_temp
    results['Min temp for Voc to reach max inverter voltage (°C)'] = temp_limit_voc_max_inv_val

    temp_limit_vmpp_upper_mppt_val = float('nan')
    if module_temp_coeff_voc != 0:
        temp_limit_vmpp_upper_mppt_val = (
            (100 * ((inverter_vmpp_max / (design_modules_per_string * module_vmpp)) - 1))
            / module_temp_coeff_voc
        ) + STC_temp
    results['Min temp for Vmpp to reach MPPT limit (upper) (°C)'] = temp_limit_vmpp_upper_mppt_val

    temp_limit_vmpp_lower_mppt_val = float('nan')
    if module_temp_coeff_voc != 0:
        temp_limit_vmpp_lower_mppt_val = (
            (100 * ((inverter_vmpp_min / (design_modules_per_string * module_vmpp)) - 1))
            / module_temp_coeff_voc
        ) + STC_temp
    results['Max temp for Vmpp to reach MPPT limit (lower) (°C)'] = temp_limit_vmpp_lower_mppt_val

    if results['array_current_change_per_degC'] != 0:
        isc_temp_difference_to_limit = results['total_current_change_max_temp'] / results['array_current_change_per_degC']
        results['Max temp for Isc to reach limit (°C)'] = isc_temp_difference_to_limit
    else:
        results['Max temp for Isc to reach limit (°C)'] = float('inf')

    results['Check: Voc Compliance'] = "OK"
    results['Comment: Voc'] = "O.K."
    if results['String Voc (V) at Max Op Temp'] > inverter_v_system_max:
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Inverter's DC input voltage is exceeded at Max Temp!"
    elif results['String Voc (V) at Min Op Temp'] > inverter_v_system_max:
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Inverter's DC input voltage is exceeded at Min Temp!"

    results['Check: Vmpp Compliance'] = "OK"
    results['Comment: Vmpp'] = "O.K."
    if results['String Vmpp (V) at Max Op Temp'] > inverter_vmpp_max:
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is exceeded at Max Temp!"
    elif results['String Vmpp (V) at Max Op Temp'] < inverter_vmpp_min:
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is under operation at Max Temp!"
    elif results['String Vmpp (V) at Min Op Temp'] < inverter_vmpp_min:
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is under operation at Min Temp!"

    results['Check: Isc Compliance'] = "OK"
    results['Comment: Isc'] = "O.K."
    if results['Array Isc (A) at Max Op Temp'] > inverter_max_pv_current_a:
        results['Check: Isc Compliance'] = "NOK"
        results['Comment: Isc'] = "Inverter's max. DC current is exceeded at Max Temp!"

    return results

# ==============================================================================
# DISPLAY FUNCTIONS (Preserved from original)
# ==============================================================================

def display_design_summary(params, results_calc):
    st.subheader("Proposed Configuration Summary")
    st.markdown(f"""
    | Parameter | Value | Unit |
    | :---------- | :---- | :--- |
    | Azimuth | {params['design_azimuth']} | degrees |
    | Tilt angle | {params['design_tilt_angle']} | degrees |
    | Row spacing | {params['design_row_spacing_m']:.2f} | m |
    | PV module rated power | {params['design_pv_module_rated_power_wp']} | Wp |
    | Ratio modules/string | {results_calc['Design Modules per String']} | |
    | Ratio strings/inverter | {results_calc['Design Strings per Inverter']} | |
    | Number of inverter | {params['design_num_inverters']} | |
    | Inverter rated AC power | {params['design_inverter_rated_ac_power_kVA']} | kVA |
    | Total rated power P$_{{DC}}$ | {results_calc['Total rated power PDC (MWp)']:.3f} | MWp |
    | Total rated power P$_{{AC}}$ | {results_calc['Total rated power PAC (MW)']:.3f} | MW |
    | Ratio P$_{{DC}}$ P$_{{AC}}$ | {results_calc['Ratio PDC PAC']:.2f} | |
    """, unsafe_allow_html=True)

def display_inverter_features(params, results_calc):
    st.subheader("Inverter ")
    st.markdown(f"""
    <div style="background-color: #e6f2ff; padding: 10px; border-radius: 5px;">
    <h5 style="margin-top: 0; color: #004d99;">Voltage Features & Current Features</h5>
    <table style="width:100%; border-collapse: collapse; text-align: center;">
      <thead>
        <tr>
          <th style="border: 1px solid #ddd; padding: 8px;"></th>
          <th style="border: 1px solid #ddd; padding: 8px;">DC V<sub>max</sub> (V) = {params['inverter_v_system_max']}</th>
          <th colspan="2" style="border: 1px solid #ddd; padding: 8px;">DC V<sub>mpp</sub> range (V) = {params['inverter_vmpp_min']} - {params['inverter_vmpp_max']}</th>
          <th colspan="2" style="border: 1px solid #ddd; padding: 8px;">DC I<sub>max</sub> (A) = {params['inverter_max_pv_current_a']}</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Physical connections</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Voc (V)</td>
          <td style="border: 1px solid #ddd; padding: 8px; background-color: #f0f8ff;">Max string</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Max string (Voc_limit_calc)']}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">1</td>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Vmpp (V)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Min string (Vmpp_min_limit_calc)']}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Max string (Vmpp_max_limit_calc)']}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Impp (A)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px; background-color: #f0f8ff;">Max string (p/MPPT)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Max string (p/MPPT)_limit_calc']}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
        </tr>
      </tbody>
    </table>
    <div style="margin-top: 10px; display: flex; justify-content: space-around; font-size: 0.9em;">
        <span style="color: green;">{results_calc['Calculated String Voc (STC) for Max String (N16)']:.1f} {results_calc['Check Voc (STC) vs Inverter Max (N17)']}</span>
        <span style="color: green;">{results_calc['Calculated String Vmpp (STC) for Min String (O16)']:.2f} {results_calc['Check Vmpp (STC) vs Inverter Min (O17)']}</span>
        <span style="color: green;">{results_calc['Calculated String Vmpp (STC) for Max String (P16)']:.2f} {results_calc['Check Vmpp (STC) vs Inverter Max (P17)']}</span>
    </div>
    </div>
    """, unsafe_allow_html=True)

def display_configuration_maximum(params, results_calc):
    st.subheader("Configuration Maximum")
    check_status = results_calc['Check_Config_DC_Power']
    color = 'green' if check_status == 'OK' else 'red'
    symbol = '✔' if check_status == 'OK' else 'X'

    st.markdown(f"""
    <div style="background-color: #e6f2ff; padding: 10px; border-radius: 5px;">
    <table style="width:100%; border-collapse: collapse; text-align: center;">
        <thead>
            <tr>
                <th style="border: 1px solid #ddd; padding: 8px;"></th>
                <th style="border: 1px solid #ddd; padding: 8px;">Number of strings per inverter</th>
                <th style="1px solid #ddd; padding: 8px;">DC Power per inverter (kWp)</th>
                <th style="border: 1px solid #ddd; padding: 8px;">Recommended DC Power per inverter (kW)</th>
                <th style="border: 1px solid #ddd; padding: 8px;">Check</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Configuration maximum</td>
                <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Configured_num_strings_per_inverter']}</td>
                <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Configured_DC_Power_per_inverter_kWp']:.2f}</td>
                <td style="border: 1px solid #ddd; padding: 8px;">{params['inverter_max_recommended_pv_power_kw']}</td>
                <td style="border: 1px solid #ddd; padding: 8px; color: {color};">{symbol}</td>
            </tr>
        </tbody>
    </table>
    </div>
    """, unsafe_allow_html=True)

def display_temperature_behaviour(params, results_calc):
    def format_value(value, decimal_places=0):
        if math.isinf(value):
            return "Inf"
        if math.isnan(value):
            return "N/A"
        return f"{value:.{decimal_places}f}"

    st.subheader("Temperature Behaviour")
    
    # Get additional temperature calculation values
    string_voc_change_per_degC = format_value(results_calc.get('string_voc_change_per_degC', float('nan')), 2)
    array_current_change_per_degC = format_value(results_calc.get('array_current_change_per_degC', float('nan')), 4)
    max_op_temp_c = format_value(params.get('max_op_temp_c', float('nan')), 0)
    
    # Get the critical temperature limits
    isc_temp_limit = format_value(results_calc.get('Max temp for Isc to reach limit (°C)', float('nan')), 0)
    
    # Additional values from your example
    total_voc_change_max_temp = format_value(results_calc.get('total_voc_change_max_temp', float('nan')), 2)
    array_impp_stc = format_value(results_calc.get('array_impp_stc', float('nan')), 2)
    array_impp_at_max_op_temp = format_value(results_calc.get('array_impp_at_max_op_temp', float('nan')), 1)
    
    st.markdown(f"""
    <div style="background-color: #e6f2ff; padding: 15px; border-radius: 5px;">
    
    <!-- Main Temperature Behaviour Table -->
    <table style="width:100%; border-collapse: collapse; text-align: center; margin-bottom: 20px;">
      <thead>
        <tr>
          <th style="border: 1px solid #ddd; padding: 8px;">Temperature Behaviour</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Temperature Coefficient</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Change per °C</th>
          <th colspan="2" style="border: 1px solid #ddd; padding: 8px;">Operating Temperature Range (°C)</th>
          <th colspan="2" style="border: 1px solid #ddd; padding: 8px;">Inverter Limits</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Check</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Comments</th>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px;"></td>
          <td style="border: 1px solid #ddd; padding: 8px;">%/°C</td>
          <td style="border: 1px solid #ddd; padding: 8px;">V/°C or A/°C</td>
          <td style="border: 1px solid #ddd; padding: 8px;">Max Temp</td>
          <td style="border: 1px solid #ddd; padding: 8px;">Min Temp</td>
          <td style="border: 1px solid #ddd; padding: 8px;">Max</td>
          <td style="border: 1px solid #ddd; padding: 8px;">Min</td>
          <td style="border: 1px solid #ddd; padding: 8px;"></td>
          <td style="border: 1px solid #ddd; padding: 8px;"></td>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Voc (V)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['module_temp_coeff_voc']:.3f}</td>
          <td style="border: 1px solid #ddd; padding: 8px; background-color: #f0f8ff;">{string_voc_change_per_degC}</td>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">{results_calc['String Voc (V) at Max Op Temp']:.0f} V</td>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold; color: green;">{results_calc['String Voc (V) at Min Op Temp']:.0f} V</td>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">{params['inverter_v_system_max']} V</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Check: Voc Compliance'] == 'OK' else 'red'};">{'✔' if results_calc['Check: Voc Compliance'] == 'OK' else 'X'}</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Comment: Voc'] == 'O.K.' else 'red'};">{results_calc['Comment: Voc']}</td>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Vmpp (V)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['module_temp_coeff_voc']:.3f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['String Vmpp (V) at Max Op Temp']:.0f} V</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['String Vmpp (V) at Min Op Temp']:.0f} V</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['inverter_vmpp_max']} V</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['inverter_vmpp_min']} V</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Check: Vmpp Compliance'] == 'OK' else 'red'};">{'✔' if results_calc['Check: Vmpp Compliance'] == 'OK' else 'X'}</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Comment: Vmpp'] == 'O.K.' else 'red'};">{results_calc['Comment: Vmpp']}</td>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Isc/Impp (A)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['module_temp_coeff_isc']:.3f}</td>
          <td style="border: 1px solid #ddd; padding: 8px; background-color: #f0f8ff;">{array_current_change_per_degC}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Array Isc (A) at Max Op Temp']:.0f} A</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Array Isc (A) at Min Op Temp']:.0f} A</td>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">{params['inverter_max_pv_current_a']} A</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Check: Isc Compliance'] == 'OK' else 'red'};">{'✔' if results_calc['Check: Isc Compliance'] == 'OK' else 'X'}</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Comment: Isc'] == 'O.K.' else 'red'};">{results_calc['Comment: Isc']}</td>
        </tr>
      </tbody>
    </table>
    
    <!-- Detailed Temperature Calculations -->
    <div style="background-color: #f8f9fa; padding: 15px; border-radius: 5px; border: 1px solid #dee2e6;">
    <h5 style="margin-top: 0; color: #004d99;">Detailed Temperature Calculations</h5>
    
    <table style="width:100%; border-collapse: collapse; text-align: left;">
      <tbody>
        <tr>
          <td style="padding: 5px; width: 40%;"><strong>String Voc change per °C:</strong></td>
          <td style="padding: 5px;">{string_voc_change_per_degC} V/°C</td>
          <td style="padding: 5px; width: 40%;"><strong>Array Current change per °C:</strong></td>
          <td style="padding: 5px;">{array_current_change_per_degC} A/°C</td>
        </tr>
        <tr>
          <td style="padding: 5px;"><strong>Max Operating Temperature:</strong></td>
          <td style="padding: 5px;">{max_op_temp_c} °C</td>
          <td style="padding: 5px;"><strong>Total Voc change at Max Temp:</strong></td>
          <td style="padding: 5px;">{total_voc_change_max_temp} V</td>
        </tr>
        <tr>
          <td style="padding: 5px;"><strong>Array Impp at STC:</strong></td>
          <td style="padding: 5px;">{array_impp_stc} A</td>
          <td style="padding: 5px;"><strong>Array Impp at Max Op Temp:</strong></td>
          <td style="padding: 5px;">{array_impp_at_max_op_temp} A</td>
        </tr>
      </tbody>
    </table>
    
    <!-- Critical Temperature Limits -->
    <div style="margin-top: 15px; display: flex; justify-content: space-between;">
      <div style="background-color: #d4edda; padding: 10px; border-radius: 5px; border: 1px solid #c3e6cb; flex: 1; margin-right: 10px;">
        <div style="font-size: 0.9em; color: #155724; margin-bottom: 3px;">Voc at Maximum Temperature</div>
        <div style="font-size: 18px; font-weight: bold; color: #155724;">{results_calc['String Voc (V) at Max Op Temp']:.0f} V</div>
      </div>
      <div style="background-color: #fff3cd; padding: 10px; border-radius: 5px; border: 1px solid #ffeaa7; flex: 1; margin-left: 10px;">
        <div style="font-size: 0.9em; color: #856404; margin-bottom: 3px;">Max temp for Isc to reach limit</div>
        <div style="font-size: 18px; font-weight: bold; color: #856404;">{isc_temp_limit} °C</div>
      </div>
    </div>
    
    </div>
    </div>
    
    <div style="margin-top: 10px; font-size: 0.9em; color: #666;">
        <strong>Note:</strong> 
        <ul style="margin-top: 5px; margin-bottom: 5px;">
            <li>String Voc change per °C: Voltage change for the entire string per degree Celsius temperature change</li>
            <li>Array Current change per °C: Current change for the entire array per degree Celsius temperature change</li>
            <li>Max temp for Isc to reach limit: Temperature at which array Isc reaches inverter's maximum current ({params['inverter_max_pv_current_a']} A)</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
       
def display_critical_temp_limits(results_calc):
    def format_value(value, decimal_places=0):
        if math.isinf(value):
            return "Inf"
        if math.isnan(value):
            return "N/A"
        return f"{value:.{decimal_places}f}"

    st.subheader("Critical Temperature Limits")
    st.markdown(f"""
    <table style="width:100%; border-collapse: collapse; text-align: left;">
        <tbody>
            <tr style="background-color: #d4edda; border: 1px solid #c3e6cb;">
                <td style="padding: 5px; font-weight: bold; color: red;">NEW</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">Min temp for Voc to reach maximum inverter voltage</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">{format_value(results_calc['Min temp for Voc to reach max inverter voltage (°C)'], 0)} °C</td>
            </tr>
            <tr style="background-color: #d4edda; border: 1px solid #c3e6cb;">
                <td style="padding: 5px; font-weight: bold; color: red;">NEW</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">Min temp for Vmpp to reach MPPT limit</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">{format_value(results_calc['Min temp for Vmpp to reach MPPT limit (upper) (°C)'], 0)} °C</td>
            </tr>
            <tr style="background-color: #d4edda; border: 1px solid #c3e6cb;">
                <td style="padding: 5px; font-weight: bold; color: red;">NEW</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">Max temp for Vmpp to reach MPPT limit</td>
                <td style="border: 1px solid #ddd; padding: 8px;">{format_value(results_calc['Max temp for Vmpp to reach MPPT limit (lower) (°C)'], 0)} °C</td>
            </tr>
        </tbody>
    </table>
    """, unsafe_allow_html=True)

# ==============================================================================
# Main Streamlit Application
# ==============================================================================

st.set_page_config(layout="wide", page_title="Solar PV Design Tool")

# Minimalist CSS
st.markdown("""
<style>
    .main-header {
        font-size: 24px;
        font-weight: 600;
        color: #333;
        margin-bottom: 20px;
        padding-bottom: 10px;
        border-bottom: 1px solid #eee;
    }
    .section-header {
        font-size: 16px;
        font-weight: 500;
        color: #444;
        margin-top: 15px;
        margin-bottom: 8px;
    }
    .compact-input {
        margin-bottom: 5px;
    }
    .stNumberInput > div > div > input {
        font-size: 14px;
        padding: 4px 8px;
    }
    .stTextInput > div > div > input {
        font-size: 14px;
        padding: 4px 8px;
    }
    .stSelectbox > div > div > select {
        font-size: 14px;
        padding: 4px 8px;
    }
    .info-text {
        font-size: 13px;
        color: #666;
        margin-top: 3px;
        margin-bottom: 8px;
    }
    .success-box {
        background-color: #f0fff4;
        border-left: 4px solid #38a169;
        padding: 10px;
        margin: 10px 0;
        border-radius: 4px;
    }
    .warning-box {
        background-color: #fffaf0;
        border-left: 4px solid #dd6b20;
        padding: 10px;
        margin: 10px 0;
        border-radius: 4px;
    }
    .temperature-display {
        background-color: #f7fafc;
        padding: 12px;
        border-radius: 6px;
        margin: 10px 0;
        border: 1px solid #e2e8f0;
    }
    .temp-value {
        font-size: 18px;
        font-weight: 600;
        color: #2d3748;
    }
    .temp-label {
        font-size: 12px;
        color: #718096;
        margin-bottom: 3px;
    }
</style>
""", unsafe_allow_html=True)

# --- Sidebar Project Manager (New) ---
with st.sidebar:
    st.markdown("### 📁 Project Library")
    
    # Save current state
    with st.expander("Save New Project"):
        p_name = st.text_input("Name", key="new_proj_name")
        if st.button("💾 Save Project"):
            if p_name:
                save_project_to_db(p_name, st.session_state.input_params)
                st.success(f"'{p_name}' saved.")
                st.rerun()
            else:
                st.error("Enter a name.")

    # Load existing project
    projects = list_projects()
    if projects:
        with st.expander("Load Project"):
            selected_p = st.selectbox("Select", projects)
            if st.button("📂 Load"):
                loaded = load_project_from_db(selected_p)
                if loaded:
                    st.session_state.input_params = loaded
                    st.session_state.show_results = False
                    st.rerun()
        
        with st.expander("Manage Database"):
            del_p = st.selectbox("Delete", projects, key="del_proj")
            if st.button("🗑️ Delete"):
                conn = sqlite3.connect(DB_NAME)
                conn.execute("DELETE FROM projects WHERE name = ?", (del_p,))
                conn.commit()
                conn.close()
                st.rerun()
    else:
        st.info("No saved projects.")

st.markdown('<div class="main-header">Solar PV Design and Validation Tool</div>', unsafe_allow_html=True)

# Initialize session state (Modified slightly to respect loaded data)
if 'input_params' not in st.session_state:
    st.session_state.input_params = {
        'module_supplier': "JA Solar", 
        'module_type': "JAM72D42-620/LB", 
        'module_vmpp': 43.51, 
        'module_voc': 51.86,
        'module_impp': 14.25, 
        'module_isc': 14.98, 
        'module_power_stc': 620, 
        'module_v_max_system': 1500,
        'module_temp_coeff_pmax': -0.300, 
        'module_temp_coeff_voc': -0.260, 
        'module_temp_coeff_isc': 0.046,
        'module_noct': 45, 
        'module_dim_width': 1134, 
        'module_dim_length': 2278,
        'selected_coeff_a': -3.47, 'selected_coeff_b': -0.0594, 'selected_coeff_delta_tcnd': 3,
        'inverter_supplier': "Sungrow", 
        'inverter_type': "SG1100UD-MV", 
        'inverter_transformer_integrated': "Yes",
        'inverter_vmpp_min': 895, 
        'inverter_vmpp_max': 1500, 
        'inverter_v_system_max': 1500,
        'inverter_max_recommended_pv_power_kw': 1650, 
        'inverter_nominal_pv_power_kw': 1100,
        'inverter_max_pv_current_a': 1435, 
        'inverter_nominal_pv_current_a': 1435,
        'inverter_nb_inputs_cc': 5, 
        'inverter_isc_max_per_inputs': 3528.0,
        'design_azimuth': 0, 'design_tilt_angle': 10, 'design_row_spacing_m': 6.00,
        'design_pv_module_rated_power_wp': 620,
        'design_modules_per_string_input': "26,28", 
        'design_strings_per_inverter_input': "90,100", 
        'design_num_inverters': 1,
        'design_inverter_rated_ac_power_kVA': 1100,
        'max_op_temp_c': float('nan'), 'min_op_temp_c': float('nan'),
        'max_temp_inclination_gain': 1.0, 'min_temp_inclination_gain': 1.0,
        'uploaded_temp_df': None,
        'last_selected_module_type': "Glass/cell/polymer sheet",
        'last_selected_mount_type': "Open rack",
        'ghi_filter_value': 50
    }

if 'show_results' not in st.session_state:
    st.session_state.show_results = False

# Main app logic - two screens
if not st.session_state.show_results:
    # INPUT SCREEN
    uploaded_file = st.file_uploader("Upload Excel file with GHI, DIF, TEMP, WS columns", type=["xlsx"])
    
    if uploaded_file is not None:
        try:
            df_temp_data = pd.read_excel(uploaded_file)
            required_cols = ['GHI', 'DIF', 'TEMP', 'WS']
            
            if all(col in df_temp_data.columns for col in required_cols):
                st.success("File uploaded successfully")
                st.session_state.input_params['uploaded_temp_df'] = df_temp_data.to_json()
                
                # Module temperature model selection
                st.markdown("### Module Temperature Model")
                
                module_types = sorted(list(set([k[0] for k in COEFFICIENT_LOOKUP.keys()])))
                mount_types = sorted(list(set([k[1] for k in COEFFICIENT_LOOKUP.keys()])))
                
                default_module_type = st.session_state.input_params.get('last_selected_module_type', "Glass/cell/polymer sheet")
                default_mount_type = st.session_state.input_params.get('last_selected_mount_type', "Open rack")
                
                col1, col2 = st.columns(2)
                with col1:
                    selected_module_type = st.selectbox("Module Type", options=module_types, 
                                                      index=module_types.index(default_module_type) if default_module_type in module_types else 0)
                with col2:
                    selected_mount_type = st.selectbox("Mount Type", options=mount_types, 
                                                      index=mount_types.index(default_mount_type) if default_mount_type in mount_types else 0)
                
                st.session_state.input_params['last_selected_module_type'] = selected_module_type
                st.session_state.input_params['last_selected_mount_type'] = selected_mount_type
                
                # Get coefficients
                lookup_key = (selected_module_type, selected_mount_type)
                if lookup_key in COEFFICIENT_LOOKUP:
                    coeffs = COEFFICIENT_LOOKUP[lookup_key]
                    st.session_state.input_params['selected_coeff_a'] = coeffs['a']
                    st.session_state.input_params['selected_coeff_b'] = coeffs['b']
                    st.session_state.input_params['selected_coeff_delta_tcnd'] = coeffs['delta_tcnd']
                    
                    col_a, col_b, col_d = st.columns(3)
                    with col_a:
                        st.write(f"**a:** {coeffs['a']:.2f}")
                    with col_b:
                        st.write(f"**b:** {coeffs['b']:.4f}")
                    with col_d:
                        st.write(f"**ΔTcnd:** {coeffs['delta_tcnd']}°C")
                
                # Inclination gain factors
                st.markdown("### Inclination Gain Factors")
                col_gain1, col_gain2 = st.columns(2)
                with col_gain1:
                    st.session_state.input_params['max_temp_inclination_gain'] = st.number_input(
                        "Max Temp Gain", 
                        value=st.session_state.input_params['max_temp_inclination_gain'], 
                        format="%.2f",
                        help="Gain factor for maximum temperature calculation"
                    )
                with col_gain2:
                    st.session_state.input_params['min_temp_inclination_gain'] = st.number_input(
                        "Min Temp Gain", 
                        value=st.session_state.input_params['min_temp_inclination_gain'], 
                        format="%.2f",
                        help="Gain factor for minimum temperature calculation"
                    )
                
                # Calculate temperatures
                a_coeff = st.session_state.input_params['selected_coeff_a']
                b_coeff = st.session_state.input_params['selected_coeff_b']
                delta_tcnd_coeff = st.session_state.input_params['selected_coeff_delta_tcnd']
                max_gain = st.session_state.input_params['max_temp_inclination_gain']
                min_gain = st.session_state.input_params['min_temp_inclination_gain']
                
                # Max Temp Calculation
                max_temp_candidates_50 = df_temp_data.sort_values(by='TEMP', ascending=False).head(50).copy()
                if not max_temp_candidates_50.empty:
                    def get_tm_value_only(row):
                        ws_adj = max(row['WS'], 0.1)
                        exp_term = EULER_E**(a_coeff + b_coeff * ws_adj)
                        return (row['GHI'] * max_gain * exp_term) + row['TEMP']
                    
                    tm_values_series = max_temp_candidates_50.apply(get_tm_value_only, axis=1)
                    if not tm_values_series.empty:
                        idx_max_tm = tm_values_series.idxmax()
                        winning_row_max = max_temp_candidates_50.loc[idx_max_tm]
                        st.session_state.input_params['max_op_temp_c'] = calculate_tcell_for_row_internal(
                            winning_row_max, a_coeff, b_coeff, delta_tcnd_coeff, max_gain
                        )
                
                # Min Temp Calculation
                overall_min_tcell_candidates = []
                ghi_filter_val = st.slider("GHI Filter Value", min_value=0, max_value=200, value=int(st.session_state.input_params.get('ghi_filter_value', 50)), 
                                           help="Filter for minimum temperature calculation")
                st.session_state.input_params['ghi_filter_value'] = ghi_filter_val
                
                min_temp_ghi_equal_df_filtered = df_temp_data[df_temp_data['GHI'] == ghi_filter_val].sort_values(by='TEMP', ascending=True)
                if not min_temp_ghi_equal_df_filtered.empty:
                    min_row_ghi_equal = min_temp_ghi_equal_df_filtered.loc[min_temp_ghi_equal_df_filtered['TEMP'].idxmin()]
                    overall_min_tcell_candidates.append(
                        calculate_tcell_for_row_internal(min_row_ghi_equal, a_coeff, b_coeff, delta_tcnd_coeff, min_gain)
                    )
                
                min_temp_ghi_geq_df_filtered = df_temp_data[df_temp_data['GHI'] >= ghi_filter_val].sort_values(by='TEMP', ascending=True)
                if not min_temp_ghi_geq_df_filtered.empty:
                    min_row_ghi_geq = min_temp_ghi_geq_df_filtered.loc[min_temp_ghi_geq_df_filtered['TEMP'].idxmin()]
                    overall_min_tcell_candidates.append(
                        calculate_tcell_for_row_internal(min_row_ghi_geq, a_coeff, b_coeff, delta_tcnd_coeff, min_gain)
                    )
                
                if overall_min_tcell_candidates:
                    st.session_state.input_params['min_op_temp_c'] = min(overall_min_tcell_candidates)
                
                # Display temperature results
                st.markdown('<div class="temperature-display">', unsafe_allow_html=True)
                col_temp1, col_temp2 = st.columns(2)
                with col_temp1:
                    st.markdown('<div class="temp-label">Max Operating Cell Temperature</div>', unsafe_allow_html=True)
                    st.markdown(f'<div class="temp-value">{format_value(st.session_state.input_params["max_op_temp_c"], 1)}°C</div>', unsafe_allow_html=True)
                with col_temp2:
                    st.markdown('<div class="temp-label">Min Operating Cell Temperature</div>', unsafe_allow_html=True)
                    st.markdown(f'<div class="temp-value">{format_value(st.session_state.input_params["min_op_temp_c"], 1)}°C</div>', unsafe_allow_html=True)
                st.markdown('</div>', unsafe_allow_html=True)
                    
            else:
                st.error(f"Missing required columns: {', '.join(required_cols)}")
                
        except Exception as e:
            st.error(f"Error reading file: {str(e)}")
    
    # COMPACT INPUT SECTIONS
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("### PV Module Specifications")
        st.session_state.input_params['module_supplier'] = st.text_input("Supplier", value=st.session_state.input_params['module_supplier'], key="mod_supplier")
        st.session_state.input_params['module_type'] = st.text_input("Type", value=st.session_state.input_params['module_type'], key="mod_type")
        st.session_state.input_params['module_vmpp'] = st.number_input("Vmpp (V)", value=st.session_state.input_params['module_vmpp'], format="%.2f", key="mod_vmpp")
        st.session_state.input_params['module_voc'] = st.number_input("Voc (V)", value=st.session_state.input_params['module_voc'], format="%.2f", key="mod_voc")
        st.session_state.input_params['module_impp'] = st.number_input("Impp (A)", value=st.session_state.input_params['module_impp'], format="%.2f", key="mod_impp")
        st.session_state.input_params['module_isc'] = st.number_input("Isc (A)", value=st.session_state.input_params['module_isc'], format="%.1f", key="mod_isc")
        st.session_state.input_params['module_power_stc'] = st.number_input("Power STC (Wp)", value=st.session_state.input_params['module_power_stc'], key="mod_power")
        st.session_state.input_params['module_v_max_system'] = st.number_input("V max (V)", value=st.session_state.input_params['module_v_max_system'], key="mod_vmax")
        st.session_state.input_params['module_temp_coeff_pmax'] = st.number_input("μPmax (%/°C)", value=st.session_state.input_params['module_temp_coeff_pmax'], format="%.3f", key="mod_tc_pmax")
        st.session_state.input_params['module_temp_coeff_voc'] = st.number_input("μVoc (%/°C)", value=st.session_state.input_params['module_temp_coeff_voc'], format="%.3f", key="mod_tc_voc")
        st.session_state.input_params['module_temp_coeff_isc'] = st.number_input("μIsc (%/°C)", value=st.session_state.input_params['module_temp_coeff_isc'], format="%.3f", key="mod_tc_isc")
        st.session_state.input_params['module_noct'] = st.number_input("NOCT (°C)", value=st.session_state.input_params['module_noct'], key="mod_noct")
    
    with col2:
        st.markdown("### Inverter Specifications")
        st.session_state.input_params['inverter_supplier'] = st.text_input("Inverter Supplier", value=st.session_state.input_params['inverter_supplier'], key="inv_supplier")
        st.session_state.input_params['inverter_type'] = st.text_input("Inverter Type", value=st.session_state.input_params['inverter_type'], key="inv_type")
        st.session_state.input_params['inverter_transformer_integrated'] = st.selectbox(
            "Transformer Integrated", 
            options=["Yes", "No"], 
            index=0 if st.session_state.input_params['inverter_transformer_integrated'] == "Yes" else 1,
            key="inv_transformer"
        )
        st.session_state.input_params['inverter_vmpp_min'] = st.number_input("Vmpp min (V)", value=st.session_state.input_params['inverter_vmpp_min'], key="inv_vmpp_min")
        st.session_state.input_params['inverter_vmpp_max'] = st.number_input("Vmpp max (V)", value=st.session_state.input_params['inverter_vmpp_max'], key="inv_vmpp_max")
        st.session_state.input_params['inverter_v_system_max'] = st.number_input("V system max (V)", value=st.session_state.input_params['inverter_v_system_max'], key="inv_vsys_max")
        st.session_state.input_params['inverter_max_recommended_pv_power_kw'] = st.number_input("Max rec. PV power (kW)", value=st.session_state.input_params['inverter_max_recommended_pv_power_kw'], key="inv_max_power")
        st.session_state.input_params['inverter_max_pv_current_a'] = st.number_input("Max PV current (A)", value=st.session_state.input_params['inverter_max_pv_current_a'], key="inv_max_current")
        st.session_state.input_params['inverter_nb_inputs_cc'] = st.number_input("Inputs CC", value=st.session_state.input_params['inverter_nb_inputs_cc'], key="inv_inputs")
        st.session_state.input_params['inverter_isc_max_per_inputs'] = st.number_input("Isc max per input", value=st.session_state.input_params['inverter_isc_max_per_inputs'], format="%.1f", key="inv_isc_max")
    
    with col3:
        st.markdown("### Design Configuration")
        st.session_state.input_params['design_azimuth'] = st.number_input("Azimuth (0°=NORTH)", value=st.session_state.input_params['design_azimuth'], key="des_azimuth")
        st.session_state.input_params['design_tilt_angle'] = st.number_input("Tilt angle (±)", value=st.session_state.input_params['design_tilt_angle'], key="des_tilt")
        st.session_state.input_params['design_row_spacing_m'] = st.number_input("Row spacing (m)", value=st.session_state.input_params['design_row_spacing_m'], format="%.2f", key="des_spacing")
        st.session_state.input_params['design_pv_module_rated_power_wp'] = st.number_input("PV module power (Wp)", value=st.session_state.input_params['design_pv_module_rated_power_wp'], key="des_power")
        st.session_state.input_params['design_modules_per_string_input'] = st.text_input(
            "Modules/string (comma-separated)",
            value=st.session_state.input_params['design_modules_per_string_input'],
            key="des_modules"
        )
        st.session_state.input_params['design_strings_per_inverter_input'] = st.text_input(
            "Strings/inverter (comma-separated)",
            value=st.session_state.input_params['design_strings_per_inverter_input'],
            key="des_strings"
        )
        st.session_state.input_params['design_num_inverters'] = st.number_input("Number of inverters", value=st.session_state.input_params['design_num_inverters'], min_value=1, key="des_num_inv")
        st.session_state.input_params['design_inverter_rated_ac_power_kVA'] = st.number_input("Inverter AC power (kVA)", value=st.session_state.input_params['design_inverter_rated_ac_power_kVA'], key="des_ac_power")

    st.markdown("---")
    calculate_col1, calculate_col2, calculate_col3 = st.columns([1, 2, 1])
    with calculate_col2:
        if st.button("Calculate Results", type="primary", use_container_width=True):
            st.session_state.show_results = True
            st.rerun()
            
else:
    # RESULTS SCREEN
    st.markdown('<div class="main-header">Design Validation Results</div>', unsafe_allow_html=True)
    
    params = st.session_state.input_params
    
    if math.isnan(params['max_op_temp_c']) or math.isnan(params['min_op_temp_c']):
        st.error("Operating temperatures are not defined. Please go back and upload climate data.")
        if st.button("Back to Input Screen"):
            st.session_state.show_results = False
            st.rerun()
        st.stop()
    
    modules_per_string_scenarios = parse_numbers_from_string(params['design_modules_per_string_input'], default_value=26)
    strings_per_inverter_scenarios = parse_numbers_from_string(params['design_strings_per_inverter_input'], default_value=90)
    
    if len(modules_per_string_scenarios) != len(strings_per_inverter_scenarios):
        st.error("Scenario count mismatch.")
        if st.button("Back to Input Screen"):
            st.session_state.show_results = False
            st.rerun()
        st.stop()
    
    all_calculated_results = []
    
    for i, (m_s, s_i) in enumerate(zip(modules_per_string_scenarios, strings_per_inverter_scenarios)):
        current_p = params.copy()
        current_p['design_modules_per_string'] = m_s
        current_p['design_strings_per_inverter'] = s_i
        
        try:
            res = calculate_solar_pv_design(
                current_p['module_supplier'], current_p['module_type'], current_p['module_vmpp'], current_p['module_voc'], 
                current_p['module_impp'], current_p['module_isc'], current_p['module_power_stc'], current_p['module_v_max_system'],
                current_p['module_temp_coeff_pmax'], current_p['module_temp_coeff_voc'], current_p['module_temp_coeff_isc'], 
                current_p['module_noct'], current_p['module_dim_width'], current_p['module_dim_length'],
                current_p['selected_coeff_a'], current_p['selected_coeff_b'], current_p['selected_coeff_delta_tcnd'],
                current_p['inverter_supplier'], current_p['inverter_type'], current_p['inverter_transformer_integrated'],
                current_p['inverter_vmpp_min'], current_p['inverter_vmpp_max'], current_p['inverter_v_system_max'],
                current_p['inverter_max_recommended_pv_power_kw'], current_p['inverter_nominal_pv_power_kw'],
                current_p['inverter_max_pv_current_a'], current_p['inverter_nominal_pv_current_a'],
                current_p['inverter_nb_inputs_cc'], current_p['inverter_isc_max_per_inputs'],
                current_p['design_azimuth'], current_p['design_tilt_angle'], current_p['design_row_spacing_m'],
                current_p['design_pv_module_rated_power_wp'], current_p['design_modules_per_string'], 
                current_p['design_strings_per_inverter'], current_p['design_num_inverters'], 
                current_p['design_inverter_rated_ac_power_kVA'], current_p['max_op_temp_c'], current_p['min_op_temp_c'],
                current_p['max_temp_inclination_gain'], current_p['min_temp_inclination_gain']
            )
            all_calculated_results.append((current_p, res))
        except Exception as e:
            all_calculated_results.append((current_p, {"Error": str(e), "Scenario Failed": True}))
    
    # Download logic
    st.subheader("Download Results")
    combined_dl = []
    for cp, rs in all_calculated_results:
        if "Scenario Failed" not in rs:
            row = {"Mod/String": cp['design_modules_per_string'], "Str/Inv": cp['design_strings_per_inverter'], **rs}
            combined_dl.append({k: ("N/A" if (isinstance(v, float) and (math.isinf(v) or math.isnan(v))) else v) for k, v in row.items()})
    
    if combined_dl:
        df_dl = pd.DataFrame(combined_dl)
        col_dl1, col_dl2 = st.columns(2)
        with col_dl1:
            st.download_button("Download CSV", df_dl.to_csv(index=False).encode('utf-8'), "results.csv", "text/csv", use_container_width=True)
        with col_dl2:
            output = io.BytesIO()
            df_dl.to_excel(output, index=False, engine='openpyxl')
            st.download_button("Download Excel", output.getvalue(), "results.xlsx", use_container_width=True)
    
    for i, (cp, rs) in enumerate(all_calculated_results):
        if "Scenario Failed" in rs: continue
        with st.expander(f"Scenario {i+1}: {cp['design_modules_per_string']}x{cp['design_strings_per_inverter']}", expanded=(i==0)):
            display_inverter_features(cp, rs)
            st.markdown("---")
            display_configuration_maximum(cp, rs)
            st.markdown("---")
            display_temperature_behaviour(cp, rs)
            st.markdown("---")
            display_critical_temp_limits(rs)

    if st.button("Back to Input Screen"):
        st.session_state.show_results = False
        st.rerun()