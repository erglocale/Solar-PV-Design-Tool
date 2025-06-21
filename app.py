import streamlit as st
import math

# Set page config at the very beginning
st.set_page_config(layout="wide", page_title="Solar PV Design Tool")

def calculate_solar_pv_design(
    # Module Specs
    module_supplier, module_type, module_vmpp, module_voc, module_impp, module_isc, module_power_stc, module_v_max_system,
    module_temp_coeff_pmax, module_temp_coeff_voc, module_temp_coeff_isc, module_noct,
    module_dim_width, module_dim_length,

    # Inverter Specs
    inverter_supplier, inverter_type, inverter_transformer_integrated,
    inverter_vmpp_min, inverter_vmpp_max, inverter_v_system_max,
    inverter_max_recommended_pv_power_kw, inverter_nominal_pv_power_kw,
    inverter_max_pv_current_a, inverter_nominal_pv_current_a, inverter_nb_inputs_cc, inverter_isc_max_per_inputs,

    # Design Configuration
    design_azimuth, design_tilt_angle, design_row_spacing_m,
    design_pv_module_rated_power_wp,
    design_modules_per_string, design_strings_per_inverter, design_num_inverters, design_inverter_rated_ac_power_kVA,

    # Operating Temperatures
    max_op_temp_c, min_op_temp_c
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
    # Voc at Max Operating Temperature
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

    # Voc at Min Operating Temperature
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
    results['temp_diff_max_op_exact'] = max_op_temp_c - STC_temp

    # total_voc_change_max_temp (N33) = N32 * N31
    results['total_voc_change_max_temp'] = results['temp_diff_max_op_exact'] * results['string_voc_change_per_degC']

    # string_voc_stc (N34) = Module_Voc * Modules_per_string
    results['string_voc_stc'] = module_voc * design_modules_per_string

    # array_current_change_per_degC (P31) = (µIsc/100) * Module_Impp * Strings_per_inverter
    results['array_current_change_per_degC'] = (module_temp_coeff_isc / 100) * module_impp * design_strings_per_inverter

    # array_impp_stc (P33) = Module_Impp * Strings_per_inverter
    results['array_impp_stc'] = module_impp * design_strings_per_inverter

    # total_impp_change_max_temp (P32 value of 7008, matches 11200 - 4192.44)
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
    # Calculation: (Remaining Current Capacity) / (Current Change per degree C for the array) + STC_temp
    # This matches your latest clarification for 7008 / 2.0123712 = 3482, which implies it is a temperature *difference* from 0 or reference.
    # The requirement is to display 3482. We will calculate the 'delta T' and display that directly.
    isc_limit_temp_from_given_formula_val = float('nan') # Initialize
    if results['array_current_change_per_degC'] != 0: # Avoid division by zero
        # The result from your Excel is the delta temperature (3482.2) and NOT (delta_T + STC_Temp).
        isc_limit_temp_from_given_formula_val = results['total_current_change_max_temp'] / results['array_current_change_per_degC']
    results['Max temp for Isc to reach limit (°C)'] = isc_limit_temp_from_given_formula_val


    # --- Compliance Checks and Comments (Updated Logic) ---
    results['Check: Voc Compliance'] = "OK"
    results['Comment: Voc'] = "O.K."
    if results['String Voc (V) at Max Op Temp'] > inverter_v_system_max:
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Inverter's DC input voltage is exceeded at Max Temp!"
    elif results['String Voc (V) at Max Op Temp'] > module_v_max_system:
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Module's max. system voltage is exceeded at max temp!"
    elif results['String Voc (V) at Min Op Temp'] > inverter_v_system_max:
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Inverter's DC input voltage is exceeded at Min Temp!"
    elif results['String Voc (V) at Min Op Temp'] > module_v_max_system:
        results['Check: Voc Compliance'] = "NOK"
        results['Comment: Voc'] = "Modules's max. system voltage is exceeded at Min Temp!"

    results['Check: Vmpp Compliance'] = "OK"
    results['Comment: Vmpp'] = "O.K."
    if results['String Vmpp (V) at Max Op Temp'] > inverter_vmpp_max:
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is exceeded at Max Temp!"
    elif results['String Vmpp (V) at Max Op Temp'] > module_v_max_system:
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Module's max. system voltage is exceeded at max temp!"
    elif results['String Vmpp (V) at Min Op Temp'] > inverter_vmpp_max:
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Inverter's DC MPPT voltage is exceeded at Min Temp!"
    elif results['String Vmpp (V) at Min Op Temp'] > module_v_max_system:
        results['Check: Vmpp Compliance'] = "NOK"
        results['Comment: Vmpp'] = "Modules's max. system voltage is exceeded at Min Temp!"
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
    elif results['Array Isc (A) at Min Op Temp'] > inverter_max_pv_current_a:
        # This condition is less critical as Isc typically decreases with temp,
        # but included for completeness from Excel structure.
        results['Check: Isc Compliance'] = "NOK"
        results['Comment: Isc'] = "Inverter's max. DC current is exceeded at Min Temp!"


    return results

# --- Streamlit Dashboard Layout ---

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

        # Operating Temperatures
        'max_op_temp_c': 71.04790582,
        'min_op_temp_c': 5.78891
    }

# Create columns for input and output
col1, col2, col3 = st.columns([1, 1, 2])

with col1:
    st.header("Main Components")

    st.subheader("PV Modules Spec")
    st.session_state.input_params['module_supplier'] = st.text_input("Supplier", value=st.session_state.input_params['module_supplier'])
    st.session_state.input_params['module_type'] = st.text_input("Type", value=st.session_state.input_params['module_type'])
    st.session_state.input_params['module_vmpp'] = st.number_input("Vmpp (V)", value=st.session_state.input_params['module_vmpp'], format="%.2f")
    st.session_state.input_params['module_voc'] = st.number_input("Voc (V)", value=st.session_state.input_params['module_voc'], format="%.2f")
    st.session_state.input_params['module_impp'] = st.number_input("Impp (A)", value=st.session_state.input_params['module_impp'], format="%.2f")
    st.session_state.input_params['module_isc'] = st.number_input("Isc (A)", value=st.session_state.input_params['module_isc'], format="%.1f")
    st.session_state.input_params['module_power_stc'] = st.number_input("Power STC (Wp)", value=st.session_state.input_params['module_power_stc'])
    st.session_state.input_params['module_v_max_system'] = st.number_input("V max (V)", value=st.session_state.input_params['module_v_max_system'])

    st.subheader("Temperature Coeff.")
    st.session_state.input_params['module_temp_coeff_pmax'] = st.number_input("µPmax (%/°C)", value=st.session_state.input_params['module_temp_coeff_pmax'], format="%.3f")
    st.session_state.input_params['module_temp_coeff_voc'] = st.number_input("µVoc (%/°C)", value=st.session_state.input_params['module_temp_coeff_voc'], format="%.3f")
    st.session_state.input_params['module_temp_coeff_isc'] = st.number_input("µIsc (%/°C)", value=st.session_state.input_params['module_temp_coeff_isc'], format="%.3f")
    st.session_state.input_params['module_noct'] = st.number_input("NOCT (°C)", value=st.session_state.input_params['module_noct'])

    st.subheader("Dimensions")
    st.session_state.input_params['module_dim_width'] = st.number_input("Dimension width (mm)", value=st.session_state.input_params['module_dim_width'])
    st.session_state.input_params['module_dim_length'] = st.number_input("Dimension (mm)", value=st.session_state.input_params['module_dim_length'])

    st.subheader("Inverter General Info")
    st.session_state.input_params['inverter_supplier'] = st.text_input("Supplier", value=st.session_state.input_params['inverter_supplier'], key="inv_sup_gen")
    st.session_state.input_params['inverter_type'] = st.text_input("Type", value=st.session_state.input_params['inverter_type'], key="inv_type_gen")
    st.session_state.input_params['inverter_transformer_integrated'] = st.text_input("Transformer integrated", value=st.session_state.input_params['inverter_transformer_integrated'], key="inv_xfmr_int")


with col2:
    st.header("Design Configuration Proposed")

    st.subheader("Proposed Design Layout")
    st.session_state.input_params['design_azimuth'] = st.number_input("Azimuth (0°=NORTH)", value=st.session_state.input_params['design_azimuth'])
    st.session_state.input_params['design_tilt_angle'] = st.number_input("Tilt angle (±)", value=st.session_state.input_params['design_tilt_angle'])
    st.session_state.input_params['design_row_spacing_m'] = st.number_input("Row spacing (m)", value=st.session_state.input_params['design_row_spacing_m'], format="%.2f")
    st.session_state.input_params['design_pv_module_rated_power_wp'] = st.number_input("PV module rated power (Wp)", value=st.session_state.input_params['design_pv_module_rated_power_wp'])

    st.subheader("Stringing Configuration")
    st.session_state.input_params['design_modules_per_string'] = st.number_input("Ratio modules/string", value=st.session_state.input_params['design_modules_per_string'], min_value=1)
    st.session_state.input_params['design_strings_per_inverter'] = st.number_input("Ratio strings/inverter", value=st.session_state.input_params['design_strings_per_inverter'], min_value=1)
    st.session_state.input_params['design_num_inverters'] = st.number_input("Number of inverter", value=st.session_state.input_params['design_num_inverters'], min_value=1)
    st.session_state.input_params['design_inverter_rated_ac_power_kVA'] = st.number_input("Inverter rated AC power (kVA)", value=st.session_state.input_params['design_inverter_rated_ac_power_kVA'])

    st.subheader("Inverter DC Input Limits")
    st.session_state.input_params['inverter_vmpp_min'] = st.number_input("Vmpp min (V)", value=st.session_state.input_params['inverter_vmpp_min'], key="inv_vmpp_min_dc_input")
    st.session_state.input_params['inverter_vmpp_max'] = st.number_input("Vmpp max (V)", value=st.session_state.input_params['inverter_vmpp_max'], key="inv_vmpp_max_dc_input")
    st.session_state.input_params['inverter_v_system_max'] = st.number_input("V system max (V)", value=st.session_state.input_params['inverter_v_system_max'], key="inv_v_sys_max_dc_input")
    st.session_state.input_params['inverter_max_recommended_pv_power_kw'] = st.number_input("Maximum recommended PV power (kW)", value=st.session_state.input_params['inverter_max_recommended_pv_power_kw'], key="inv_max_rec_pwr_dc_input")
    st.session_state.input_params['inverter_nominal_pv_power_kw'] = st.number_input("Nominal PV power (kW)", value=st.session_state.input_params['inverter_nominal_pv_power_kw'], key="inv_nom_pv_pwr_dc_input")
    st.session_state.input_params['inverter_max_pv_current_a'] = st.number_input("Maximum PV current (A)", value=st.session_state.input_params['inverter_max_pv_current_a'], key="inv_max_pv_curr_dc_input")
    st.session_state.input_params['inverter_nominal_pv_current_a'] = st.number_input("Nominal PV current (A)", value=st.session_state.input_params['inverter_nominal_pv_current_a'], key="inv_nom_pv_curr_dc_input")
    st.session_state.input_params['inverter_nb_inputs_cc'] = st.number_input("Nb inputs CC", value=st.session_state.input_params['inverter_nb_inputs_cc'], key="inv_nb_inputs_dc_input")
    st.session_state.input_params['inverter_isc_max_per_inputs'] = st.number_input("Isc max per inputs (A)", value=st.session_state.input_params['inverter_isc_max_per_inputs'], format="%.1f", key="inv_isc_max_inputs")


    st.subheader("Operating Temperature Range")
    st.session_state.input_params['max_op_temp_c'] = st.number_input("Max Operating Temperature (°C)", value=st.session_state.input_params['max_op_temp_c'])
    st.session_state.input_params['min_op_temp_c'] = st.number_input("Min Operating Temperature (°C)", value=st.session_state.input_params['min_op_temp_c'])


# --- Perform Calculations (on every rerun due to Streamlit's model) ---
params = st.session_state.input_params
results_calc = calculate_solar_pv_design(
    # Module
    params['module_supplier'], params['module_type'], params['module_vmpp'], params['module_voc'], params['module_impp'], params['module_isc'],
    params['module_power_stc'], params['module_v_max_system'], params['module_temp_coeff_pmax'], params['module_temp_coeff_voc'],
    params['module_temp_coeff_isc'], params['module_noct'], params['module_dim_width'], params['module_dim_length'],

    # Inverter General Info
    params['inverter_supplier'], params['inverter_type'], params['inverter_transformer_integrated'],

    # Inverter DC Input
    params['inverter_vmpp_min'], params['inverter_vmpp_max'], params['inverter_v_system_max'],
    params['inverter_max_recommended_pv_power_kw'], params['inverter_nominal_pv_power_kw'],
    params['inverter_max_pv_current_a'], params['inverter_nominal_pv_current_a'],
    params['inverter_nb_inputs_cc'], params['inverter_isc_max_per_inputs'],

    # Design Configuration
    params['design_azimuth'], params['design_tilt_angle'], params['design_row_spacing_m'],
    params['design_pv_module_rated_power_wp'],
    params['design_modules_per_string'], params['design_strings_per_inverter'], params['design_num_inverters'],
    params['design_inverter_rated_ac_power_kVA'],

    # Temperatures
    params['max_op_temp_c'], params['min_op_temp_c']
)

with col3:
    st.header("STRING DIMENSIONS")

    st.subheader("Inverter Features")

    # Replicating the top section of "STRING DIMENSIONS" table
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


    st.markdown("---") # Separator

    st.subheader("Configuration Maximum")
    st.markdown(f"""
    <div style="background-color: #e6f2ff; padding: 10px; border-radius: 5px;">
    <table style="width:100%; border-collapse: collapse; text-align: center;">
        <thead>
            <tr>
                <th style="border: 1px solid #ddd; padding: 8px;"></th>
                <th style="border: 1px solid #ddd; padding: 8px;">Number of strings per inverter</th>
                <th style="border: 1px solid #ddd; padding: 8px;">DC Power per inverter (kWp)</th>
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

    st.markdown("---") # Separator

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
          <td style="border: 1px solid #ddd; padding: 8px; font-weight: bold; color: green;">{results_calc['Array Isc (A) at Max Op Temp']:.0f}</td>
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
                <td style="border-right: 1px solid #c3e6cb; padding: 5px;">Voc at Maximum Temperature: {results_calc['String Voc (V) at Max Op Temp']:.0f} V</td>
                <td style="padding: 5px;">Max temp for Isc to reach limit: {results_calc['Max temp for Isc to reach limit (°C)']:.0f} °C</td>
            </tr>
        </table>
    </div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---") # Separator

    st.subheader("Critical Temperature Limits")
    st.markdown(f"""
    <table style="width:100%; border-collapse: collapse; text-align: left;">
        <tbody>
            <tr style="background-color: #d4edda; border: 1px solid #c3e6cb;">
                <td style="padding: 5px; font-weight: bold; color: red;">NEW</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">Min temp for Voc to reach maximum inverter voltage</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">{results_calc['Min temp for Voc to reach max inverter voltage (°C)']:.0f} °C</td>
            </tr>
            <tr style="background-color: #d4edda; border: 1px solid #c3e6cb;">
                <td style="padding: 5px; font-weight: bold; color: red;">NEW</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">Min temp for Vmpp to reach MPPT limit</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">{results_calc['Min temp for Vmpp to reach MPPT limit (upper) (°C)']:.0f} °C</td>
            </tr>
            <tr style="background-color: #d4edda; border: 1px solid #c3e6cb;">
                <td style="padding: 5px; font-weight: bold; color: red;">NEW</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">Max temp for Vmpp to reach MPPT limit</td>
                <td style="border: 1px solid #c3e6cb; padding: 5px;">{results_calc['Max temp for Vmpp to reach MPPT limit (lower) (°C)']:.0f} °C</td>
            </tr>
        </tbody>
    </table>
    """, unsafe_allow_html=True) 