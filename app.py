import streamlit as st
import math
import pandas as pd
import io

# Define Euler's Constant
EULER_E = 2.718281828

# --- Define Lookup Table for a, b, DeltaTcnd (GLOBAL SCOPE) ---
COEFFICIENT_LOOKUP = {
    ("Glass/cell/glass", "Open rack"): {"a": -3.47, "b": -0.0594, "delta_tcnd": 3},
    ("Glass/cell/glass", "Close-roof mount"): {"a": -2.98, "b": -0.0471, "delta_tcnd": 1},
    ("Glass/cell/polymer sheet", "Open rack"): {"a": -3.56, "b": -0.0750, "delta_tcnd": 3},
    ("Glass/cell/polymer sheet", "Insulated back"): {"a": -2.81, "b": -0.0455, "delta_tcnd": 0},
    ("Polymer/thin-film/steel", "Open rack"): {"a": -3.58, "b": -0.1130, "delta_tcnd": 3},
}

# --- Helper function for conditional formatting of values (GLOBAL SCOPE) ---
def format_value(value, decimal_places=0):
    if math.isinf(value):
        return "Inf"  # Return "Inf" for infinity
    if math.isnan(value):
        return "N/A"  # Return "N/A" for NaN
    return f"{value:.{decimal_places}f}"

# --- Helper function to calculate Tm and Tcell for a single row based on a given gain (GLOBAL SCOPE) ---
def calculate_tcell_for_row_internal(row, a, b, delta_tcnd, gain):
    # Ensure WS is not too low for exp function if b is negative and large
    ws_val = row['WS']
    ws_adj = max(ws_val, 0.1) if ws_val is not None else 0.1 # Ensure non-zero/non-negative for exp calculation

    # Handle GHI as well, if it's zero or None, the solar heating term should be zero
    ghi_val = row['GHI'] if row['GHI'] is not None else 0

    # Prevent math.exp from causing overflow/underflow if a+b*ws is extreme (very large positive/negative)
    exp_term_arg = a + b * ws_adj
    if exp_term_arg > 700: # Approx ln(float_max) to prevent overflow
        exp_term = float('inf')
    elif exp_term_arg < -700: # Approx ln(float_min) to prevent underflow
        exp_term = 0.0
    else:
        exp_term = EULER_E**(exp_term_arg)

    # Tm = (GHI * Inclination Gain * (e^(a + (b * WS)))) + Temp
    tm_val = (ghi_val * gain * exp_term) + row['TEMP']
    # Tcell = Tm + ((GHI * Inclination Gain)/1000) * delta_tcnd
    tcell_val = tm_val + ((ghi_val * gain) / 1000) * delta_tcnd

    return tcell_val


# Set page config at the very beginning
st.set_page_config(layout="wide", page_title="Solar PV Design Tool")

# --- Core Calculation Function ---
# Renamed module_v_system to module_v_max_system and inverter_v_system to inverter_v_system_max to match input_params keys
def calculate_solar_pv_design(
    # Module Specs
    module_supplier, module_type, module_vmpp, module_voc, module_impp, module_isc, module_power_stc, module_v_max_system,
    module_temp_coeff_pmax, module_temp_coeff_voc, module_temp_coeff_isc, module_noct,
    module_dim_width, module_dim_length,

    # Added a, b, delta_tcnd from dropdown selection
    selected_coeff_a, selected_coeff_b, selected_coeff_delta_tcnd,

    # Inverter Specs
    inverter_supplier, inverter_type, inverter_transformer_integrated,
    inverter_vmpp_min, inverter_vmpp_max, inverter_v_system_max, # Ensure this is inverter_v_system_max consistently
    inverter_max_recommended_pv_power_kw, inverter_nominal_pv_power_kw,
    inverter_max_pv_current_a, inverter_nominal_pv_current_a, inverter_nb_inputs_cc, inverter_isc_max_per_inputs,

    # Design Configuration
    design_azimuth, design_tilt_angle, design_row_spacing_m,
    design_pv_module_rated_power_wp,
    design_modules_per_string, design_strings_per_inverter, design_num_inverters, design_inverter_rated_ac_power_kVA,

    # Operating Temperatures - Derived from uploaded data in Step 4
    max_op_temp_c, min_op_temp_c, # These are the *final derived Tcell values* from Step 4
    max_temp_inclination_gain, min_temp_inclination_gain # Inclination gains are inputs
):
    """
    Performs all solar PV string design and validation calculations.

    Args:
        (All parameters directly from Streamlit inputs)

    Returns:
        dict: A dictionary containing all calculated results, checks, and comments.
    """

    results = {}
    STC_temp = 25  # Standard Test Condition temperature

    # Initialize temperature limit keys to prevent KeyError if calculation conditions are not met
    results['Min temp for Voc to reach max inverter voltage (°C)'] = float('nan')
    results['Min temp for Vmpp to reach MPPT limit (upper) (°C)'] = float('nan')
    results['Max temp for Vmpp to reach MPPT limit (lower) (°C)'] = float('nan')
    results['Max temp for Isc to reach limit (°C)'] = float('nan')


    # --- Design Configuration Proposed Calculations ---
    total_rated_power_PDC_Wp = (
        design_modules_per_string * design_strings_per_inverter * design_num_inverters * design_pv_module_rated_power_wp
    )
    results['Total rated power PDC (MWp)'] = total_rated_power_PDC_Wp / 1_000_000

    total_rated_power_PAC_MW = design_inverter_rated_ac_power_kVA / 1000
    results['Total rated power PAC (MW)'] = total_rated_power_PAC_MW

    if total_rated_power_PAC_MW != 0:
        results['Ratio PDC PAC'] = results['Total rated power PDC (MWp)'] / total_rated_power_PAC_MW
    else:
        results['Ratio PDC PAC'] = 0 # Handle division by zero

    # --- STRING DIMENSIONS - TOP SECTION CALCULATIONS ---

    # N11: Max string (Voc limit)
    theoretical_max_modules_voc = inverter_v_system_max / module_voc
    rounded_max_modules_voc = round(theoretical_max_modules_voc)
    if rounded_max_modules_voc * module_voc > inverter_v_system_max:
        max_modules_voc_calc = rounded_max_modules_voc - 1
    else:
        max_modules_voc_calc = rounded_max_modules_voc
    results['Max string (Voc_limit_calc)'] = max_modules_voc_calc # M11

    # O12: Min string (Vmpp_min_limit)
    theoretical_min_modules_vmpp = inverter_vmpp_min / module_vmpp
    rounded_min_modules_vmpp = round(theoretical_min_modules_vmpp)
    if rounded_min_modules_vmpp * module_vmpp < inverter_vmpp_min:
        min_modules_vmpp_calc = rounded_min_modules_vmpp + 1
    else:
        min_modules_vmpp_calc = rounded_min_modules_vmpp
    results['Min string (Vmpp_min_limit_calc)'] = min_modules_vmpp_calc # O12

    # P12: Max string (Vmpp_max_limit)
    theoretical_max_modules_vmpp = inverter_vmpp_max / module_vmpp
    rounded_max_modules_vmpp = round(theoretical_max_modules_vmpp)
    if rounded_max_modules_vmpp * module_vmpp > inverter_vmpp_max:
        max_modules_vmpp_calc = rounded_max_modules_vmpp - 1
    else:
        max_modules_vmpp_calc = rounded_max_modules_vmpp
    results['Max string (Vmpp_max_limit_calc)'] = max_modules_vmpp_calc # P12

    # R13: Max string (p/MPPT) from I_max_inv / Impp_module
    if module_impp != 0:
        max_strings_per_mppt_limit_calc = math.floor(inverter_max_pv_current_a / module_impp)
    else:
        max_strings_per_mppt_limit_calc = 0
    results['Max string (p/MPPT)_limit_calc'] = max_strings_per_mppt_limit_calc # R13

    # Intermediate String Voltage/Current Checks (Below the Max String row in Excel)
    results['Calculated String Voc (STC) for Max String (N16)'] = max_modules_voc_calc * module_voc
    results['Check Voc (STC) vs Inverter Max (N17)'] = "OK" if results['Calculated String Voc (STC) for Max String (N16)'] < inverter_v_system_max else "NOK"

    results['Calculated String Vmpp (STC) for Min String (O16)'] = min_modules_vmpp_calc * module_vmpp
    results['Check Vmpp (STC) vs Inverter Min (O17)'] = "OK" if results['Calculated String Vmpp (STC) for Min String (O16)'] > inverter_vmpp_min else "NOK"

    results['Calculated String Vmpp (STC) for Max String (P16)'] = max_modules_vmpp_calc * module_vmpp
    results['Check Vmpp (STC) vs Inverter Max (P17)'] = "OK" if results['Calculated String Vmpp (STC) for Max String (P16)'] < inverter_vmpp_max else "NOK"


    # --- Configuration Maximum Table ---
    results['Configured_num_strings_per_inverter'] = design_strings_per_inverter # N21
    dc_power_per_inverter_kWp = (design_modules_per_string * design_strings_per_inverter * module_power_stc) / 1000
    results['Configured_DC_Power_per_inverter_kWp'] = dc_power_per_inverter_kWp # O21

    results['Check_Config_DC_Power'] = "OK" if dc_power_per_inverter_kWp <= inverter_max_recommended_pv_power_kw else "NOK" # S21


    # --- Temperature Behavior Calculations ---
    # Voc at Max Operating Temperature (using the max_op_temp_c provided, which is max T_cell from Step 4)
    temp_diff_max_op = max_op_temp_c - STC_temp
    voc_temp_factor_max = 1 + (module_temp_coeff_voc * temp_diff_max_op / 100)
    string_voc_at_max_op_temp = design_modules_per_string * module_voc * voc_temp_factor_max
    results['String Voc (V) at Max Op Temp'] = string_voc_at_max_op_temp

    # Vmpp at Max Operating Temperature (using Voc coeff as per Excel sheet)
    vmpp_temp_factor_max = 1 + (module_temp_coeff_voc * temp_diff_max_op / 100)
    string_vmpp_at_max_op_temp = design_modules_per_string * module_vmpp * vmpp_temp_factor_max
    results['String Vmpp (V) at Max Op Temp'] = string_vmpp_at_max_op_temp

    # Isc at Max Operating Temperature (Array Total)
    isc_temp_factor_max = 1 + (module_temp_coeff_isc * temp_diff_max_op / 100)
    array_isc_at_max_op_temp = (
        design_strings_per_inverter * module_isc * isc_temp_factor_max
    )
    results['Array Isc (A) at Max Op Temp'] = array_isc_at_max_op_temp

    # Voc at Min Operating Temperature (using the min_op_temp_c provided, which is min T_cell from Step 4)
    temp_diff_min_op = min_op_temp_c - STC_temp
    voc_temp_factor_min = 1 + (module_temp_coeff_voc * temp_diff_min_op / 100)
    string_voc_at_min_op_temp = design_modules_per_string * module_voc * voc_temp_factor_min
    results['String Voc (V) at Min Op Temp'] = string_voc_at_min_op_temp

    # Vmpp at Min Operating Temperature (using Voc coeff as per Excel sheet)
    vmpp_temp_factor_min = 1 + (module_temp_coeff_voc * temp_diff_min_op / 100)
    string_vmpp_at_min_op_temp = design_modules_per_string * module_vmpp * vmpp_temp_factor_min
    results['String Vmpp (V) at Min Op Temp'] = string_vmpp_at_min_op_temp

    # Isc at Min Operating Temperature (Array Total)
    isc_temp_factor_min = 1 + (module_temp_coeff_isc * temp_diff_min_op / 100)
    array_isc_at_min_op_temp = (
        design_strings_per_inverter * module_isc * isc_temp_factor_min
    )
    results['Array Isc (A) at Min Op Temp'] = array_isc_at_min_op_temp


    # --- Additional Green Outputs Calculations (From the image) ---
    # string_voc_change_per_degC (N31) = (µVoc/100) * Module_Voc * Modules_per_string
    results['string_voc_change_per_degC'] = (module_temp_coeff_voc / 100) * module_voc * design_modules_per_string

    # temp_diff_max_op_exact (N32) = Max_Op_Temp - 25
    results['temp_diff_max_op_exact'] = max_op_temp_c - STC_temp # Uses the derived max_op_temp_c

    # total_voc_change_max_temp (N33) = N32 * N31
    results['total_voc_change_max_temp'] = results['temp_diff_max_op_exact'] * results['string_voc_change_per_degC']

    # string_voc_stc (N34) = Module_Voc * Modules_per_string
    results['string_voc_stc'] = module_voc * design_modules_per_string

    # array_current_change_per_degC (P31) = (µIsc/100) * Module_Impp * Strings_per_inverter
    results['array_current_change_per_degC'] = (module_temp_coeff_isc / 100) * module_impp * design_strings_per_inverter

    # array_impp_stc (P33) = Module_Impp * Strings_per_inverter
    results['array_impp_stc'] = module_impp * design_strings_per_inverter

    # total_impp_change_max_temp (P32 value of 7008, matches inverter_max_pv_current_a - array_impp_stc)
    results['total_current_change_max_temp'] = inverter_max_pv_current_a - results['array_impp_stc']

    # array_impp_at_max_op_temp (P34 in Excel, likely current at max temp)
    results['array_impp_at_max_op_temp'] = results['array_impp_stc'] + (results['temp_diff_max_op_exact'] * results['array_current_change_per_degC'])


    # --- Temperature Limits for Compliance (Q37, Q38, Q39, and P35) ---
    # Min temp for Voc to reach maximum inverter voltage (Q37)
    temp_limit_voc_max_inv_val = float('nan') # Initialize
    if module_temp_coeff_voc != 0:
        temp_limit_voc_max_inv_val = (
            (100 * ((inverter_v_system_max / (design_modules_per_string * module_voc)) - 1))
            / module_temp_coeff_voc
        ) + STC_temp
    results['Min temp for Voc to reach max inverter voltage (°C)'] = temp_limit_voc_max_inv_val


    # Min temp for Vmpp to reach MPPT limit (upper bound) (Q38)
    temp_limit_vmpp_upper_mppt_val = float('nan') # Initialize
    if module_temp_coeff_voc != 0: # Using Voc coeff for Vmpp calculation, as in Excel.
        temp_limit_vmpp_upper_mppt_val = (
            (100 * ((inverter_vmpp_max / (design_modules_per_string * module_vmpp)) - 1))
            / module_temp_coeff_voc # If Vmpp coeff is separate, use that.
        ) + STC_temp
    results['Min temp for Vmpp to reach MPPT limit (upper) (°C)'] = temp_limit_vmpp_upper_mppt_val


    # Max temp for Vmpp to reach MPPT limit (lower bound) (Q39)
    temp_limit_vmpp_lower_mppt_val = float('nan') # Initialize
    if module_temp_coeff_voc != 0: # Using Voc coeff for Vmpp calculation, as in Excel.
        temp_limit_vmpp_lower_mppt_val = (
            (100 * ((inverter_vmpp_min / (design_modules_per_string * module_vmpp)) - 1))
            / module_temp_coeff_voc # If Vmpp coeff is separate, use that.
        ) + STC_temp
    results['Max temp for Vmpp to reach MPPT limit (lower) (°C)'] = temp_limit_vmpp_lower_mppt_val

    # Max temp for Isc to reach limit (P35 in Excel context)
    isc_limit_temp_from_given_formula_val = float('nan') # Initialize
    if results['array_current_change_per_degC'] != 0: # Avoid division by zero
        isc_temp_difference_to_limit = results['total_current_change_max_temp'] / results['array_current_change_per_degC']
        results['Max temp for Isc to reach limit (°C)'] = isc_temp_difference_to_limit
    else:
        results['Max temp for Isc to reach limit (°C)'] = float('inf')


    # --- Compliance Checks and Comments (Updated Logic) ---
    results['Check: Voc Compliance'] = "OK"
    results['Comment: Voc'] = "O.K."
    if results['String Voc (V) at Max Op Temp'] > inverter_v_system_max: # Direct argument
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Inverter's DC input voltage is exceeded at Max Temp!"
    elif results['String Voc (V) at Max Op Temp'] > module_v_max_system: # Direct argument
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Module's max. system voltage is exceeded at max temp!"
    elif results['String Voc (V) at Min Op Temp'] > inverter_v_system_max: # Direct argument
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Inverter's DC input voltage is exceeded at Min Temp!"
    elif results['String Voc (V) at Min Op Temp'] > module_v_max_system: # Direct argument
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Modules's max. system voltage is exceeded at Min Temp!"

    results['Check: Vmpp Compliance'] = "OK"
    results['Comment: Vmpp'] = "O.K."
    if results['String Vmpp (V) at Max Op Temp'] > inverter_vmpp_max: # Direct argument
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is exceeded at Max Temp!"
    elif results['String Vmpp (V) at Max Op Temp'] > module_v_max_system: # Direct argument
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Module's max. system voltage is exceeded at max temp!"
    elif results['String Vmpp (V) at Min Op Temp'] > inverter_vmpp_max: # Direct argument
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is exceeded at Min Temp!"
    elif results['String Vmpp (V) at Min Op Temp'] > module_v_max_system: # Direct argument
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Modules's max. system voltage is exceeded at Min Temp!"
    elif results['String Vmpp (V) at Max Op Temp'] < inverter_vmpp_min: # Direct argument
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is under operation at Max Temp!"
    elif results['String Vmpp (V) at Min Op Temp'] < inverter_vmpp_min: # Direct argument
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is under operation at Min Temp!"


    results['Check: Isc Compliance'] = "OK"
    results['Comment: Isc'] = "O.K."
    if results['Array Isc (A) at Max Op Temp'] > inverter_max_pv_current_a: # Direct argument
        results['Check: Isc Compliance'] = "NOK"
        results['Comment: Isc'] = "Inverter's max. DC current is exceeded at Max Temp!"
    elif results['Array Isc (A) at Min Op Temp'] > inverter_max_pv_current_a: # Direct argument
        # This condition is less critical as Isc typically decreases with temp,
        # but included for completeness from Excel structure.
        results['Check: Isc Compliance'] = "NOK"
        results['Comment: Isc'] = "Inverter's max. DC current is exceeded at Min Temp!"


    return results

# --- Streamlit Display Functions (Refactored) ---

def display_design_summary(params, results_calc):
    st.subheader("Proposed Configuration Summary")
    st.markdown(f"""
    | Parameter | Value | Unit |
    | :---------- | :---- | :--- |
    | Azimuth | {params['design_azimuth']} | degrees |
    | Tilt angle | {params['design_tilt_angle']} | degrees |
    | Row spacing | {params['design_row_spacing_m']:.2f} | m |
    | PV module rated power | {params['design_pv_module_rated_power_wp']} | Wp |
    | Ratio modules/string | {params['design_modules_per_string']} | |
    | Ratio strings/inverter | {params['design_strings_per_inverter']} | |
    | Number of inverter | {params['design_num_inverters']} | |
    | Inverter rated AC power | {params['design_inverter_rated_ac_power_kVA']} | kVA |
    | Total rated power P$_{{DC}}$ | {results_calc['Total rated power PDC (MWp)']:.3f} | MWp |
    | Total rated power P$_{{AC}}$ | {results_calc['Total rated power PAC (MW)']:.3f} | MW |
    | Ratio P$_{{DC}}$ P$_{{AC}}$ | {results_calc['Ratio PDC PAC']:.2f} | |
    """, unsafe_allow_html=True)

def display_inverter_features(params, results_calc):
    st.subheader("Inverter Features")
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
    <div style="margin-top: 5px; font-size: 0.9em; text-align: center;">
        <p>String size: [{results_calc['Min string (Vmpp_min_limit_calc)']}-{results_calc['Max string (Voc_limit_calc)']}] mod/string;</p>
        <p>Maximum inverter load: {results_calc['Max string (p/MPPT)_limit_calc']} strings</p>
    </div>
    </div>
    """, unsafe_allow_html=True)

def display_configuration_maximum(params, results_calc):
    st.subheader("Configuration Maximum")
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
                <td style="border: 1px solid #ddd; padding: 8px; color: green;">✔</td>
            </tr>
        </tbody>
    </table>
    </div>
    """, unsafe_allow_html=True)

def display_temperature_behaviour(params, results_calc):
    # Helper function for conditional formatting of values (NaN/inf to "N/A")
    def format_value(value, decimal_places=0):
        if math.isinf(value):
            return "Inf"  # Return "Inf" for infinity
        if math.isnan(value):
            return "N/A"  # Return "N/A" for NaN
        return f"{value:.{decimal_places}f}"

    st.subheader("Temperature Behaviour")

    st.markdown(f"""
    <div style="background-color: #e6f2ff; padding: 10px; border-radius: 5px;">
    <table style="width:100%; border-collapse: collapse; text-align: center;">
      <thead>
        <tr>
          <th style="border: 1px solid #ddd; padding: 8px;">Temperature Behaviour</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Temperature Coefficient</th>
          <th colspan="2" style="border: 1px solid #ddd; padding: 8px;">Operating Temperature Range (°C)</th>
          <th colspan="2" style="border: 1px solid #ddd; padding: 8px;">Inverter range</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Check</th>
          <th style="border: 1px solid #ddd; padding: 8px;">Comments</th>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px;"></td>
          <td style="border: 1px solid #ddd; padding: 8px;">%/°C</td>
          <td style="border: 1px solid #ddd; padding: 8px;">Max</td>
          <td style="border: 1px solid #ddd; padding: 8px;">Min</td>
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
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['String Voc (V) at Max Op Temp']:.0f}</td>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold; color: green;">{results_calc['String Voc (V) at Min Op Temp']:.0f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['inverter_v_system_max']}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Check: Voc Compliance'] == 'OK' else 'red'};">{'✔' if results_calc['Check: Voc Compliance'] == 'OK' else 'X'}</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Comment: Voc'] == 'O.K.' else 'red'};">{results_calc['Comment: Voc']}</td>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Vmpp (V)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['module_temp_coeff_voc']:.3f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['String Vmpp (V) at Max Op Temp']:.0f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['String Vmpp (V) at Min Op Temp']:.0f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['inverter_vmpp_max']}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['inverter_vmpp_min']}</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Check: Vmpp Compliance'] == 'OK' else 'red'};">{'✔' if results_calc['Check: Vmpp Compliance'] == 'OK' else 'X'}</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Comment: Vmpp'] == 'O.K.' else 'red'};">{results_calc['Comment: Vmpp']}</td>
        </tr>
        <tr>
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold;">Isc (A)</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['module_temp_coeff_isc']:.3f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Array Isc (A) at Max Op Temp']:.0f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{results_calc['Array Isc (A) at Min Op Temp']:.0f}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">{params['inverter_max_pv_current_a']}</td>
          <td style="border: 1px solid #ddd; padding: 8px;">-</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Check: Isc Compliance'] == 'OK' else 'red'};">{'✔' if results_calc['Check: Isc Compliance'] == 'OK' else 'X'}</td>
          <td style="border: 1px solid #ddd; padding: 8px; color: {'green' if results_calc['Comment: Isc'] == 'O.K.' else 'red'};">{results_calc['Comment: Isc']}</td>
        </tr>
      </tbody>
    </table>
    <div style="margin-top: 10px; display: flex; justify-content: space-around; font-size: 0.9em;">
        <span style="color: green;">{results_calc['string_voc_change_per_degC']:.5f} V/°C</span>
        <span style="color: green;">{results_calc['array_current_change_per_degC']:.7f} A/°C</span>
    </div>
    <div style="display: flex; justify-content: space-around; font-size: 0.9em;">
        <span style="color: green;">{results_calc['temp_diff_max_op_exact']:.0f}</span>
        <span style="color: green;">{results_calc['total_current_change_max_temp']:.2f}</span>
    </div>
    <div style="display: flex; justify-content: space-around; font-size: 0.9em;">
        <span style="color: green;">{results_calc['total_voc_change_max_temp']:.5f}</span>
        <span style="color: green;">{results_calc['array_impp_stc']:.2f}</span>
    </div>
    <div style="display: flex; justify-content: space-around; font-size: 0.9em;">
        <span style="color: green;">{results_calc['string_voc_stc']:.1f}</span>
        <span style="color: green;">{params['inverter_max_pv_current_a']}</span>
    </div>
    <div style="margin-top: 15px; background-color: #d4edda; border: 1px solid #c3e6cb; padding: 10px; border-radius: 5px;">
        <table style="width:100%; border-collapse: collapse; text-align: left;">
            <tr>
                <td style="border-right: 1px solid #c3e6cb; padding: 5px;">Voc at Maximum Temperature: {format_value(results_calc['String Voc (V) at Max Op Temp'], 0)} V</td>
                <td style="padding: 5px;">Max temp for Isc to reach limit: {format_value(results_calc['Max temp for Isc to reach limit (°C)'], 0)} °C</td>
            </tr>
        </table>
    </div>
    </div>
    """, unsafe_allow_html=True)

def display_critical_temp_limits(results_calc):
    # Helper function for conditional formatting of values (NaN/inf to "N/A")
    def format_value(value, decimal_places=0):
        if math.isinf(value):
            return "Inf"  # Return "Inf" for infinity
        if math.isnan(value):
            return "N/A"  # Return "N/A" for NaN
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
                <td style="border: 1px solid #c3e6cb; padding: 5px;">{format_value(results_calc['Max temp for Vmpp to reach MPPT limit (lower) (°C)'], 0)} °C</td>
            </tr>
        </tbody>
    </table>
    """, unsafe_allow_html=True)

# --- Main Streamlit App Execution Flow ---
st.title("Solar PV Design and Validation Tool")
st.markdown("Use this tool to design your PV array stringing and check for electrical compatibility.")

# Initialize session state for inputs with all defaults from the images
if 'input_params' not in st.session_state:
    st.session_state.input_params = {
        # PV Modules Spec
        'module_supplier': "Jinko",
        'module_type': "530",
        'module_vmpp': 40.71,
        'module_voc': 49.35,
        'module_impp': 13.02,
        'module_isc': 13.7,
        'module_power_stc': 530,
        'module_v_max_system': 1500,
        'module_temp_coeff_pmax': -0.350,
        'module_temp_coeff_voc': -0.280,
        'module_temp_coeff_isc': 0.048,
        'module_noct': 42,
        'module_dim_width': 1134,
        'module_dim_length': 2274,

        # Default a, b, delta_tcnd for initial run if no selection made
        'selected_coeff_a': -3.47,
        'selected_coeff_b': -0.0594,
        'selected_coeff_delta_tcnd': 3,


        # Inverter General Info
        'inverter_supplier': "SMA",
        'inverter_type': "Sunny Cnetral 4400",
        'inverter_transformer_integrated': "N/A", # Not specified in image

        # Inverter DC input
        'inverter_vmpp_min': 962,
        'inverter_vmpp_max': 1325,
        'inverter_v_system_max': 1500,
        'inverter_max_recommended_pv_power_kw': 7000,
        'inverter_nominal_pv_power_kw': 55000,
        'inverter_max_pv_current_a': 11200,
        'inverter_nominal_pv_current_a': 8400,
        'inverter_nb_inputs_cc': 32,
        'inverter_isc_max_per_inputs': 0.0, # Not specified in image

        # Design Configuration Proposed Inputs
        'design_azimuth': 0,
        'design_tilt_angle': 60,
        'design_row_spacing_m': 6.35,
        'design_pv_module_rated_power_wp': 530,
        'design_modules_per_string': 28,
        'design_strings_per_inverter': 322,
        'design_num_inverters': 1,
        'design_inverter_rated_ac_power_kVA': 3960,

        # Operating Temperatures - Will be set by upload process or default NaN
        'max_op_temp_c': float('nan'),
        'min_op_temp_c': float('nan'),
        'max_temp_inclination_gain': 1.0, # New input
        'min_temp_inclination_gain': 1.3, # New input
        # New session state variable to store the uploaded DataFrame itself
        'uploaded_temp_df': None
    }

# --- Step-by-step navigation ---
# Total steps: 0 (Welcome), 1 (Module), 2 (Inverter), 3 (Design), 4 (Temp Range), 5 (Results)
# So, max step index is 5.
if 'current_step' not in st.session_state:
    st.session_state.current_step = 0

# Navigation buttons
col_nav1, col_nav2, col_nav3 = st.columns([1, 1, 1])

with col_nav1:
    if st.session_state.current_step > 0:
        if st.button("⬅️ Previous Step", use_container_width=True):
            st.session_state.current_step -= 1
            st.rerun()

with col_nav3: # Place Next button on the right
    if st.session_state.current_step < 5: # If not on the final results step
        if st.button("Next Step ➡️", use_container_width=True):
            st.session_state.current_step += 1
            st.rerun()
    elif st.session_state.current_step == 5: # On the results step, offer to recalculate
        if st.button("Recalculate / View Results Again", use_container_width=True):
            st.session_state.current_step = 0 # Go back to start or a specific step
            st.rerun()

st.markdown("---") # Separator below navigation buttons


# --- Conditional Step Content ---
params = st.session_state.input_params

# Define columns for input steps where applicable
# These will be instantiated inside each relevant if block.
# This avoids NameError by instantiating columns ONLY when needed for a given step.
# And ensures content inside these columns are correctly aligned.

if st.session_state.current_step == 0:
    st.header("Welcome!")
    st.write("Click 'Next Step' to begin configuring your solar PV system.")
    st.write("You will go through different input sections before viewing the final results.")

elif st.session_state.current_step == 1:
    main_col_left, main_col_right = st.columns([1, 1]) # Columns specific to this step
    with main_col_left:
        st.header("Step 1: PV Module Specifications")
        st.subheader("PV Modules Spec")
        st.session_state.input_params['module_supplier'] = st.text_input("Supplier", value=params['module_supplier'])
        st.session_state.input_params['module_type'] = st.text_input("Type", value=params['module_type'])
        st.session_state.input_params['module_vmpp'] = st.number_input("Vmpp (V)", value=params['module_vmpp'], format="%.2f")
        st.session_state.input_params['module_voc'] = st.number_input("Voc (V)", value=params['module_voc'], format="%.2f")
        st.session_state.input_params['module_impp'] = st.number_input("Impp (A)", value=params['module_impp'], format="%.2f")
        st.session_state.input_params['module_isc'] = st.number_input("Isc (A)", value=params['module_isc'], format="%.1f")
        st.session_state.input_params['module_power_stc'] = st.number_input("Power STC (Wp)", value=params['module_power_stc'])
        st.session_state.input_params['module_v_max_system'] = st.number_input("V max (V)", value=params['module_v_max_system'])

    with main_col_right: # Place temperature coeffs and dimensions in right column
        st.subheader("Temperature Coeff.")
        st.session_state.input_params['module_temp_coeff_pmax'] = st.number_input("µPmax (%/°C)", value=params['module_temp_coeff_pmax'], format="%.3f")
        st.session_state.input_params['module_temp_coeff_voc'] = st.number_input("µVoc (%/°C)", value=params['module_temp_coeff_voc'], format="%.3f")
        st.session_state.input_params['module_temp_coeff_isc'] = st.number_input("µIsc (%/°C)", value=params['module_temp_coeff_isc'], format="%.3f")
        st.session_state.input_params['module_noct'] = st.number_input("NOCT (°C)", value=params['module_noct'])

        st.subheader("Dimensions")
        st.session_state.input_params['module_dim_width'] = st.number_input("Dimension width (mm)", value=params['module_dim_width'])
        st.session_state.input_params['module_dim_length'] = st.number_input("Dimension (mm)", value=params['module_dim_length'])


elif st.session_state.current_step == 2:
    main_col_left, main_col_right = st.columns([1, 1]) # Columns specific to this step
    with main_col_left:
        st.header("Step 2: Inverter Specifications")
        st.subheader("Inverter General Info")
        st.session_state.input_params['inverter_supplier'] = st.text_input("Supplier", value=params['inverter_supplier'])
        st.session_state.input_params['inverter_type'] = st.text_input("Type", value=params['inverter_type'])
        st.session_state.input_params['inverter_transformer_integrated'] = st.text_input("Transformer integrated", value=params['inverter_transformer_integrated'])

    with main_col_right:
        st.subheader("Inverter DC Input Limits")
        st.session_state.input_params['inverter_vmpp_min'] = st.number_input("Vmpp min (V)", value=params['inverter_vmpp_min'], key="inv_vmpp_min_dc_input")
        st.session_state.input_params['inverter_vmpp_max'] = st.number_input("Vmpp max (V)", value=params['inverter_vmpp_max'], key="inv_vmpp_max_dc_input")
        st.session_state.input_params['inverter_v_system_max'] = st.number_input("V system max (V)", value=params['inverter_v_system_max'], key="inv_v_sys_max_dc_input")
        st.session_state.input_params['inverter_max_recommended_pv_power_kw'] = st.number_input("Maximum recommended PV power (kW)", value=params['inverter_max_recommended_pv_power_kw'], key="inv_max_rec_pwr_dc_input")
        st.session_state.input_params['inverter_nominal_pv_power_kw'] = st.number_input("Nominal PV power (kW)", value=params['inverter_nominal_pv_power_kw'], key="inv_nom_pv_pwr_dc_input")
    # Continue with current inputs in right column
        st.session_state.input_params['inverter_max_pv_current_a'] = st.number_input("Maximum PV current (A)", value=params['inverter_max_pv_current_a'], key="inv_max_pv_curr_dc_input")
        st.session_state.input_params['inverter_nominal_pv_current_a'] = st.number_input("Nominal PV current (A)", value=params['inverter_nominal_pv_current_a'], key="inv_nom_pv_curr_dc_input")
        st.session_state.input_params['inverter_nb_inputs_cc'] = st.number_input("Nb inputs CC", value=params['inverter_nb_inputs_cc'], key="inv_nb_inputs_dc_input")
        st.session_state.input_params['inverter_isc_max_per_inputs'] = st.number_input("Isc max per inputs (A)", value=params['inverter_isc_max_per_inputs'], format="%.1f", key="inv_isc_max_inputs")


elif st.session_state.current_step == 3:
    main_col_left, main_col_right = st.columns([1, 1]) # Columns specific to this step
    with main_col_left:
        st.header("Step 3: Design Configuration")
        st.subheader("Proposed Design Layout")
        st.session_state.input_params['design_azimuth'] = st.number_input("Azimuth (0°=NORTH)", value=params['design_azimuth'])
        st.session_state.input_params['design_tilt_angle'] = st.number_input("Tilt angle (±)", value=params['design_tilt_angle'])
        st.session_state.input_params['design_row_spacing_m'] = st.number_input("Row spacing (m)", value=params['design_row_spacing_m'], format="%.2f")
        st.session_state.input_params['design_pv_module_rated_power_wp'] = st.number_input("PV module rated power (Wp)", value=params['design_pv_module_rated_power_wp'])

    with main_col_right:
        st.subheader("Stringing Configuration")
        st.session_state.input_params['design_modules_per_string'] = st.number_input("Ratio modules/string", value=params['design_modules_per_string'], min_value=1)
        st.session_state.input_params['design_strings_per_inverter'] = st.number_input("Ratio strings/inverter", value=params['design_strings_per_inverter'], min_value=1)
        st.session_state.input_params['design_num_inverters'] = st.number_input("Number of inverter", value=params['design_num_inverters'], min_value=1)
        st.session_state.input_params['design_inverter_rated_ac_power_kVA'] = st.number_input("Inverter rated AC power (kVA)", value=params['design_inverter_rated_ac_power_kVA'])


elif st.session_state.current_step == 4: # Dedicated step for Operating Temperature Range
    st.header("Step 4: Operating Temperature Range")

    # --- Module Temperature Model (a, b, ΔTcnd) ---
    st.subheader("Module Temperature Model (Coefficients)")
    st.write("These coefficients are used to predict cell temperature: Tcell = Tm + ((GHI * Inclination Gain)/1000) * ΔTcnd, where Tm = (GHI * Inclination Gain * e^(a + b * WS)) + Temp")

    # Dropdowns for Module Type and Mount
    module_types = sorted(list(set([k[0] for k in COEFFICIENT_LOOKUP.keys()]))) # Get unique module types
    mount_types = sorted(list(set([k[1] for k in COEFFICIENT_LOOKUP.keys()]))) # Get unique mount types

    # Find the specific default indices for "Glass/cell/polymer sheet" and "Open rack"
    default_module_type_str = "Glass/cell/polymer sheet"
    default_mount_type_str = "Open rack"
    default_module_idx = module_types.index(default_module_type_str) if default_module_type_str in module_types else 0
    default_mount_idx = mount_types.index(default_mount_type_str) if default_mount_type_str in mount_types else 0


    # Use current selection if it exists, otherwise use the specific default index
    current_mod_type_idx = module_types.index(params.get('last_selected_module_type', default_module_type_str)) if params.get('last_selected_module_type') in module_types else default_module_idx
    current_mount_type_idx = mount_types.index(params.get('last_selected_mount_type', default_mount_type_str)) if params.get('last_selected_mount_type') in mount_types else default_mount_idx


    col_model_type, col_model_mount = st.columns(2)
    with col_model_type:
        selected_module_type = st.selectbox("Select Module Type", options=module_types, index=current_mod_type_idx, key="sel_mod_type")
    with col_model_mount:
        selected_mount_type = st.selectbox("Select Mount Type", options=mount_types, index=current_mount_type_idx, key="sel_mount_type")

    # Store last selections for persistence across reruns
    st.session_state.input_params['last_selected_module_type'] = selected_module_type
    st.session_state.input_params['last_selected_mount_type'] = selected_mount_type


    # Lookup and display a, b, DeltaTcnd
    lookup_key = (selected_module_type, selected_mount_type)
    if lookup_key in COEFFICIENT_LOOKUP:
        coeffs = COEFFICIENT_LOOKUP[lookup_key]
        st.session_state.input_params['selected_coeff_a'] = coeffs['a']
        st.session_state.input_params['selected_coeff_b'] = coeffs['b']
        st.session_state.input_params['selected_coeff_delta_tcnd'] = coeffs['delta_tcnd']
        st.write(f"**a:** {coeffs['a']:.2f}")
        st.write(f"**b:** {coeffs['b']:.4f}")
        st.write(f"**ΔTcnd (°C):** {coeffs['delta_tcnd']:.0f}")
    else:
        st.warning("Coefficients for this Module Type and Mount combination not found. Defaulting to (Glass/cell/polymer sheet, Open rack).")
        # Default to the specific default coeffs if selected combination isn't found
        default_coeffs = COEFFICIENT_LOOKUP[(default_module_type_str, default_mount_type_str)]
        st.session_state.input_params['selected_coeff_a'] = default_coeffs['a']
        st.session_state.input_params['selected_coeff_b'] = default_coeffs['b']
        st.session_state.input_params['selected_coeff_delta_tcnd'] = default_coeffs['delta_tcnd']
        st.write(f"**a (default):** {st.session_state.input_params['selected_coeff_a']:.2f}")
        st.write(f"**b (default):** {st.session_state.input_params['selected_coeff_b']:.4f}")
        st.write(f"**ΔTcnd (°C) (default):** {st.session_state.input_params['selected_coeff_delta_tcnd']:.0f}")

    st.markdown("---")
    st.subheader("Inclination Gain Factors")
    col_gain1, col_gain2 = st.columns(2)
    with col_gain1:
        st.session_state.input_params['max_temp_inclination_gain'] = st.number_input("Max Temp Inclination Gain", value=params['max_temp_inclination_gain'], format="%.1f", key="max_temp_gain")
    with col_gain2:
        st.session_state.input_params['min_temp_inclination_gain'] = st.number_input("Min Temp Inclination Gain", value=params['min_temp_inclination_gain'], format="%.1f", key="min_temp_gain")

    st.markdown("---")
    st.subheader("Upload Environmental Data")

    # Use a unique key for the file uploader for session state
    uploaded_file = st.file_uploader("Upload Excel file (.xlsx) with GHI, DIF, TEMP, WS data", type=["xlsx"], key="temp_data_uploader")

    # Initialize df_processing as an empty DataFrame to prevent NameError if no file is valid
    df_processing = pd.DataFrame()
    valid_file_uploaded_flag = False # Flag to indicate if a valid DF is available


    # Try to load existing df from session state if available and valid
    if st.session_state.input_params['uploaded_temp_df'] is not None:
        try:
            df_temp_data_from_session = pd.read_json(st.session_state.input_params['uploaded_temp_df'])
            if not df_temp_data_from_session.empty:
                df_processing = df_temp_data_from_session.copy() # Use this as the base if valid
                valid_file_uploaded_flag = True
            else:
                st.session_state.input_params['uploaded_temp_df'] = None # Clear if it was empty
        except Exception: # Catch if JSON parsing fails (e.g., corrupted state)
            st.session_state.input_params['uploaded_temp_df'] = None # Clear bad state
            df_temp_data_from_session = None

    # Process newly uploaded file, if any
    if uploaded_file is not None:
        # Check if new file or if session state file is invalid/different
        if not valid_file_uploaded_flag or uploaded_file.name != st.session_state.get('last_uploaded_filename_temp_data'):
            try:
                df_temp_data = pd.read_excel(uploaded_file)
                st.session_state.last_uploaded_filename_temp_data = uploaded_file.name # Store filename for next rerun

                required_cols = ['GHI', 'DIF', 'TEMP', 'WS']
                if not all(col in df_temp_data.columns for col in required_cols):
                    st.error(f"Error: Missing required columns. Please ensure the file has columns: {', '.join(required_cols)}")
                    st.session_state.input_params['uploaded_temp_df'] = None
                    valid_file_uploaded_flag = False
                else:
                    df_temp_data['GHI'] = pd.to_numeric(df_temp_data['GHI'], errors='coerce')
                    df_temp_data['DIF'] = pd.to_numeric(df_temp_data['DIF'], errors='coerce')
                    df_temp_data['TEMP'] = pd.to_numeric(df_temp_data['TEMP'], errors='coerce')
                    df_temp_data['WS'] = pd.to_numeric(df_temp_data['WS'], errors='coerce')

                    if df_temp_data[['GHI', 'DIF', 'TEMP', 'WS']].isnull().any().any():
                        st.error("Error: Non-numeric data found in GHI, DIF, TEMP, or WS columns. Please ensure these columns contain only numbers.")
                        st.session_state.input_params['uploaded_temp_df'] = None
                        valid_file_uploaded_flag = False
                    else:
                        st.success("File uploaded and validated successfully!")
                        st.session_state.input_params['uploaded_temp_df'] = df_temp_data.to_json() # Store DF as JSON
                        df_processing = df_temp_data.copy() # Assign the newly uploaded and validated DF
                        valid_file_uploaded_flag = True
            except Exception as e:
                st.error(f"An error occurred while processing the file: {e}")
                st.info("Please ensure it's a valid .xlsx Excel file with correct columns and numeric data.")
                st.session_state.input_params['uploaded_temp_df'] = None
                valid_file_uploaded_flag = False
        # else: uploaded_file is not None and it matches session state, valid_file_uploaded_flag is already True
            # df_processing is already loaded from session state at the beginning of this block

    # If no valid file is uploaded yet (first load or user cleared file)
    if not valid_file_uploaded_flag:
        st.info("Please upload an Excel file to derive operating temperatures.")
        # Ensure temps are reset if no file is validly selected
        st.session_state.input_params['max_op_temp_c'] = float('nan')
        st.session_state.input_params['min_op_temp_c'] = float('nan')


    # --- Calculate and Store T_cell values for max_op_temp_c and min_op_temp_c ---
    if valid_file_uploaded_flag and not df_processing.empty:
        a_coeff = params['selected_coeff_a']
        b_coeff = params['selected_coeff_b']
        delta_tcnd_coeff = params['selected_coeff_delta_tcnd']
        max_gain = params['max_temp_inclination_gain']
        min_gain = params['min_temp_inclination_gain']

        # Max Temp Calculation: Find the max Tcell from the top 20 rows sorted by Ambient TEMP (descending)
        max_temp_df_ambient_sorted_20 = df_processing.sort_values(by='TEMP', ascending=False).head(20).copy()

        # Calculate Tcell for these 20 rows using MAX GAIN
        tcell_values_for_max_selection = max_temp_df_ambient_sorted_20.apply(
            lambda row: calculate_tcell_for_row_internal(row, a_coeff, b_coeff, delta_tcnd_coeff, max_gain),
            axis=1
        )

        if not tcell_values_for_max_selection.empty:
            max_op_temp_c_calculated = tcell_values_for_max_selection.max()
            st.session_state.input_params['max_op_temp_c'] = max_op_temp_c_calculated
        else:
            st.session_state.input_params['max_op_temp_c'] = float('nan')


        # Min Temp Calculation: Find the min Tcell from specific filtered rows using MIN GAIN
        overall_min_tcell_candidates = []

        # Retrieve the current GHI filter value (from display input or default)
        ghi_filter_value = st.session_state.get('ghi_filter_value', 50)

        # Candidate 1: Tcell for min ambient temp from GHI == user_input_GHI filter
        min_temp_ghi_equal_df_filtered = df_processing[df_processing['GHI'] == ghi_filter_value].sort_values(by='TEMP', ascending=True)
        if not min_temp_ghi_equal_df_filtered.empty:
            min_row_ghi_equal = min_temp_ghi_equal_df_filtered.loc[min_temp_ghi_equal_df_filtered['TEMP'].idxmin()]
            overall_min_tcell_candidates.append(
                calculate_tcell_for_row_internal(min_row_ghi_equal, a_coeff, b_coeff, delta_tcnd_coeff, min_gain)
            )

        # Candidate 2: Tcell for min ambient temp from GHI >= user_input_GHI filter
        min_temp_ghi_geq_df_filtered = df_processing[df_processing['GHI'] >= ghi_filter_value].sort_values(by='TEMP', ascending=True)
        if not min_temp_ghi_geq_df_filtered.empty:
            min_row_ghi_geq = min_temp_ghi_geq_df_filtered.loc[min_temp_ghi_geq_df_filtered['TEMP'].idxmin()]
            overall_min_tcell_candidates.append(
                calculate_tcell_for_row_internal(min_row_ghi_geq, a_coeff, b_coeff, delta_tcnd_coeff, min_gain)
            )

        # Select the absolute minimum among candidates
        if overall_min_tcell_candidates:
            st.session_state.input_params['min_op_temp_c'] = min(overall_min_tcell_candidates)
        else:
            st.session_state.input_params['min_op_temp_c'] = float('nan')


    # Display derived temps (even if NaN, use format_value to show N/A)
    st.info(f"Predicted Max Operating Cell Temp: **{format_value(st.session_state.input_params['max_op_temp_c'], 2)} °C**")
    st.info(f"Predicted Min Operating Cell Temp: **{format_value(st.session_state.input_params['min_op_temp_c'], 2)} °C**")


    # --- Displaying Filtered Tables and Extracted Rows (ONLY if df_processing has data) ---
    if valid_file_uploaded_flag and not df_processing.empty: # Check if df_processing has data
        st.markdown("---")
        st.subheader("Temperature Data Analysis Tables")

        # Input for GHI filter (re-defined here to be local to this block)
        ghi_filter_value = st.number_input("Enter GHI value for Min Temp filtering (e.g., 50)", min_value=0, value=st.session_state.get('ghi_filter_value', 50), step=1, key="ghi_filter_input_display")
        st.session_state.ghi_filter_value = ghi_filter_value # Store for persistence

        # Max Temp: Extract first 20 rows after sorting by TEMP (descending)
        st.markdown("#### Max Temperature Data (Top 20 by Ambient Temp)")
        max_temp_df_sorted_20 = df_processing.sort_values(by='TEMP', ascending=False).head(20)
        st.dataframe(max_temp_df_sorted_20[['Date', 'Time', 'GHI', 'DIF', 'TEMP', 'WS']]) # Display original columns only

        # Extract the row for the overall max Tcell (from the specific top 20 dataframe)
        if not max_temp_df_ambient_sorted_20.empty:
            # Calculate Tcell for all rows in this specific 20-row df with max_gain
            tcell_all_rows_for_max_gain_series = max_temp_df_ambient_sorted_20.apply(
                lambda row: calculate_tcell_for_row_internal(row, a_coeff, b_coeff, delta_tcnd_coeff, max_gain),
                axis=1
            )
            # Find the row corresponding to the maximum Tcell within *this specific 20-row dataframe*
            idx_max_tcell_in_top20 = tcell_values_for_max_selection.idxmax()
            max_tcell_overall_row_from_data = max_temp_df_ambient_sorted_20.loc[idx_max_tcell_in_top20]
            max_tcell_val_for_display = tcell_values_for_max_selection.max() # CORRECTED LINE

            st.markdown("##### Row with Max Predicted Cell Temperature from Top 20 Ambient Data")
            st.dataframe(max_tcell_overall_row_from_data[['Date', 'Time', 'GHI', 'DIF', 'TEMP', 'WS']].to_frame().T)
            st.write(f"**Predicted Cell Temp for this row:** {max_tcell_val_for_display:.2f} °C")
        else:
            st.info("No data available for Max Temperature analysis.")


        st.markdown("#### Min Temperature Data Filters")

        # Min Temp: GHI == user_input_GHI
        st.markdown(f"##### Min Temperature Data (GHI == {ghi_filter_value} W/m²)")
        min_temp_ghi_equal_df = df_processing[df_processing['GHI'] == ghi_filter_value].sort_values(by='TEMP', ascending=True)
        if not min_temp_ghi_equal_df.empty:
            st.dataframe(min_temp_ghi_equal_df[['Date', 'Time', 'GHI', 'DIF', 'TEMP', 'WS']])
            min_temp_ghi_equal_row = min_temp_ghi_equal_df.loc[min_temp_ghi_equal_df['TEMP'].idxmin()]
            st.markdown(f"###### Row with Min Ambient Temp (GHI == {ghi_filter_value})")
            st.dataframe(min_temp_ghi_equal_row[['Date', 'Time', 'GHI', 'DIF', 'TEMP', 'WS']].to_frame().T)
            st.write(f"**Predicted Cell Temp for this row:** {calculate_tcell_for_row_internal(min_temp_ghi_equal_row, a_coeff, b_coeff, delta_tcnd_coeff, min_gain):.2f} °C")
        else:
            st.info(f"No data where GHI is exactly {ghi_filter_value} W/m².")

        # Min Temp: GHI >= user_input_GHI
        st.markdown(f"##### Min Temperature Data (GHI >= {ghi_filter_value} W/m²)")
        min_temp_ghi_geq_df = df_processing[df_processing['GHI'] >= ghi_filter_value].sort_values(by='TEMP', ascending=True)
        if not min_temp_ghi_geq_df.empty:
            st.dataframe(min_temp_ghi_geq_df[['Date', 'Time', 'GHI', 'DIF', 'TEMP', 'WS']])
            min_temp_ghi_geq_row = min_temp_ghi_geq_df.loc[min_temp_ghi_geq_df['TEMP'].idxmin()]
            st.markdown(f"###### Row with Min Ambient Temp (GHI >= {ghi_filter_value})")
            st.dataframe(min_temp_ghi_geq_row[['Date', 'Time', 'GHI', 'DIF', 'TEMP', 'WS']].to_frame().T)
            st.write(f"**Predicted Cell Temp for this row:** {calculate_tcell_for_row_internal(min_temp_ghi_geq_row, a_coeff, b_coeff, delta_tcnd_coeff, min_gain):.2f} °C")
        else:
            st.info(f"No data where GHI is >= {ghi_filter_value} W/m².")

        # Removed the "Row with Overall Min Ambient Temperature" display based on explicit request.


elif st.session_state.current_step == 5: # Results step (now step 5)
    st.header("Step 5: Design Validation Results")

    # Check if temperatures were derived successfully (i.e., file was uploaded and valid)
    if math.isnan(params['max_op_temp_c']) or math.isnan(params['min_op_temp_c']):
        st.error("Operating temperatures are not yet defined. Please go back to 'Step 4: Operating Temperature Range' and upload a valid Excel file.")
        st.stop() # Stop execution if critical inputs are missing

    # --- Perform Calculations (Only on the results step, now that we know temps are set) ---
    params_for_calc = st.session_state.input_params
    results_calc = calculate_solar_pv_design(
        # Module
        params_for_calc['module_supplier'], params_for_calc['module_type'], params_for_calc['module_vmpp'], params_for_calc['module_voc'], params_for_calc['module_impp'], params_for_calc['module_isc'],
        params_for_calc['module_power_stc'], params_for_calc['module_v_max_system'], params_for_calc['module_temp_coeff_pmax'], params_for_calc['module_temp_coeff_voc'],
        params_for_calc['module_temp_coeff_isc'], params_for_calc['module_noct'], params_for_calc['module_dim_width'], params_for_calc['module_dim_length'],

        # a, b, delta_tcnd (these are passed through from Step 4)
        params_for_calc['selected_coeff_a'], params_for_calc['selected_coeff_b'], params_for_calc['selected_coeff_delta_tcnd'],

        # Inverter General Info
        params_for_calc['inverter_supplier'], params_for_calc['inverter_type'], params_for_calc['inverter_transformer_integrated'],

        # Inverter DC Input
        params_for_calc['inverter_vmpp_min'], params_for_calc['inverter_vmpp_max'], params_for_calc['inverter_v_system_max'],
        params_for_calc['inverter_max_recommended_pv_power_kw'], params_for_calc['inverter_nominal_pv_power_kw'],
        params_for_calc['inverter_max_pv_current_a'], params_for_calc['inverter_nominal_pv_current_a'],
        params_for_calc['inverter_nb_inputs_cc'], params_for_calc['inverter_isc_max_per_inputs'],

        # Design Configuration
        params_for_calc['design_azimuth'], params_for_calc['design_tilt_angle'], params_for_calc['design_row_spacing_m'],
        params_for_calc['design_pv_module_rated_power_wp'],
        params_for_calc['design_modules_per_string'], params_for_calc['design_strings_per_inverter'], params_for_calc['design_num_inverters'],
        params_for_calc['design_inverter_rated_ac_power_kVA'],

        # Temperatures (Derived from uploaded file in Step 4 + Inclination Gain)
        params_for_calc['max_op_temp_c'], params_for_calc['min_op_temp_c'],
        params_for_calc['max_temp_inclination_gain'], params_for_calc['min_temp_inclination_gain']
    )

    # --- Download Buttons ---
    st.subheader("Download Data")
    col_dl1, col_dl2, col_dl3, col_dl4 = st.columns(4)

    # Prepare Input Data for download
    input_df = pd.DataFrame(st.session_state.input_params.items(), columns=['Parameter', 'Value'])
    csv_input = input_df.to_csv(index=False).encode('utf-8')
    excel_input_buffer = io.BytesIO()
    input_df.to_excel(excel_input_buffer, index=False, engine='openpyxl')
    excel_input_buffer.seek(0)

    # Prepare Output Data for download
    display_results_for_download = {k: ("N/A" if (isinstance(v, float) and (math.isinf(v) or math.isnan(v))) else v) for k, v in results_calc.items()}
    results_df = pd.DataFrame(display_results_for_download.items(), columns=['Parameter', 'Value'])
    csv_results = results_df.to_csv(index=False).encode('utf-8')
    excel_results_buffer = io.BytesIO()
    results_df.to_excel(excel_results_buffer, index=False, engine='openpyxl')
    excel_results_buffer.seek(0)

    with col_dl1:
        st.download_button(
            label="Download Input Data (CSV)",
            data=csv_input,
            file_name="solar_input_data.csv",
            mime="text/csv",
            key="download_input_csv"
        )
    with col_dl2:
        st.download_button(
            label="Download Input Data (Excel)",
            data=excel_input_buffer,
            file_name="solar_input_data.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="download_input_excel"
        )
    with col_dl3:
        st.download_button(
            label="Download Results (CSV)",
            data=csv_results,
            file_name="solar_results.csv",
            mime="text/csv",
            key="download_results_csv"
        )
    with col_dl4:
        st.download_button(
            label="Download Results (Excel)",
            data=excel_results_buffer,
            file_name="solar_results.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="download_results_excel"
        )
    st.markdown("---") # Separator


    # --- Display Results ---
    display_design_summary(params_for_calc, results_calc)
    st.markdown("---") # Separator

    display_inverter_features(params_for_calc, results_calc)
    st.markdown("---") # Separator

    display_configuration_maximum(params_for_calc, results_calc)
    st.markdown("---") # Separator

    display_temperature_behaviour(params_for_calc, results_calc)
    st.markdown("---") # Separator

    display_critical_temp_limits(results_calc)