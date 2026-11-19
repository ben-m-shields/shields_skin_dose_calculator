def read_RDSR_irradiation_events(filename):
    import rdsr_navigator as nav
    import pandas as pd
    from pydicom import dcmread
    rdsr_obj = nav.read_file(filename)
    dicom = dcmread(filename, force=True)
    df = pd.DataFrame({
    "KVP":pd.Series([],dtype="float64"),
    "Exposure":pd.Series([],dtype="float64"),
    "Pulse_Width":pd.Series([],dtype="float64"),
    "Focal_Spot_Size":pd.Series([],dtype="float64"),
    "Pulse_Rate":pd.Series([],dtype="float64"),
    "X-Ray_Filter_Thickness_Mimimum":pd.Series([],dtype="float64"),
    "Irradiation_Event_Type":pd.Series([],dtype="string"),
    "X-Ray_Filter_Material":pd.Series([],dtype="string"),
    "Acquisition_Protocol":pd.Series([],dtype="string"),
    "Acquisition_Plane_in_Irradiation_Event":pd.Series([],dtype="string"),
    "Irradiation_Event_UID":pd.Series([],dtype="string"),
    "Positioner_Primary_Angle":pd.Series([],dtype="float64"),
    "Positioner_Secondary_Angle":pd.Series([],dtype="float64"),
    "Table_Lateral_Position":pd.Series([],dtype="float64"),
    "Table_Longitudinal_Position":pd.Series([],dtype="float64"),
    "Table_Height_Position":pd.Series([],dtype="float64"),
    "Dose_Area_Product":pd.Series([],dtype="float64"),
    "Dose_RP":pd.Series([],dtype="float64"),
    "Distance_Source_to_Isocenter":pd.Series([],dtype="float64"),
    "Distance_Source_to_Detector":pd.Series([],dtype="float64"),
    "Positioner_Primary_End_Angle":pd.Series([],dtype="float64"),
    "Positioner_Secondary_End_Angle":pd.Series([],dtype="float64")})
    
    for index, irr_event in enumerate(rdsr_obj.get_all('irradiation_event_x-ray_data')):
        try:
          df.at[index, "KVP"] = irr_event['kvp'].value[0]
        except:
          df.at[index, "KVP"] = float("nan")
        try:
          df.at[index, "Exposure"] = irr_event['exposure'].value[0]
        except:
          df.at[index, "Exposure"] = float("nan")
        try:
          df.at[index, "Pulse_Width"] = irr_event['pulse_width'].value[0]
        except:
          df.at[index, "Pulse_Width"] = float("nan")
        try:
          df.at[index, "Number of Pulses"] = irr_event['number_of_pulses'].value[0]
        except:
          df.at[index, "Number of Pulses"] = float("nan")
        try:
          df.at[index, "Focal_Spot_Size"] = irr_event['focal_spot_size'].value[0]
        except:
          df.at[index, "Focal_Spot_Size"] = float("nan")
        try:
          df.at[index, "Pulse_Rate"] = irr_event['pulse_rate'].value[0]
        except:
          df.at[index, "Pulse_Rate"] = float("nan")
        try:
          df.at[index, "Irradiation_Event_Type"] = irr_event['irradiation_event_type'].value
        except:
          df.at[index, "Irradiation_Event_Type"] = "Not provided in RDSR"
        try:
          df.at[index, "Acquisition_Protocol"] = irr_event['acquisition_protocol'].value
        except:
          df.at[index, "Acquisition_Protocol"] = "Not provided in RDSR"
        df.at[index, "Acquisition_Plane_in_Irradiation_Event"] = irr_event['acquisition_plane'].value
        df.at[index, "Irradiation_Event_UID"] = irr_event['irradiation_event_uid'].value
        df.at[index, "Positioner_Primary_Angle"] = irr_event['positioner_primary_angle'].value[0]
        df.at[index, "Positioner_Secondary_Angle"] = irr_event['positioner_secondary_angle'].value[0]
        df.at[index, "Table_Lateral_Position"] = irr_event['table_lateral_position'].value[0]
        df.at[index, "Table_Longitudinal_Position"] = irr_event['table_longitudinal_position'].value[0]
        df.at[index, "Table_Height_Position"] = irr_event['table_height_position'].value[0]
        df.at[index, "Dose_Area_Product"] = irr_event['dose_area_product'].value[0]
        df.at[index, "Dose_RP"] = irr_event['dose_(rp)'].value[0]
        try:
          df.at[index, "Distance_Source_to_Isocenter"] = irr_event['distance_source_to_isocenter'].value[0]
        except:
          df.at[index, "Distance_Source_to_Isocenter"] = float("nan")
        try:
          df.at[index, "Distance_Source_to_Detector"] = irr_event['distance_source_to_detector'].value[0]
        except:
          df.at[index, "Distance_Source_to_Detector"] = float("nan")
        try:
            df.at[index, "X-Ray_Filter_Thickness_Minimum"] = irr_event['x-ray_filters', 'x-ray_filter_thickness_minimum'].value[0]
            df.at[index, "X-Ray_Filter_Material"] = irr_event['x-ray_filters', 'x-ray_filter_material'].value
        except:
            df.at[index, "X-Ray_Filter_Thickness_Minimum"] = float("nan")
            df.at[index, "X-Ray_Filter_Material"] = "Not provided in RDSR"
        try:
            df.at[index, "Positioner_Primary_End_Angle"] = irr_event['positioner_primary_end_angle'].value[0]
            df.at[index, "Positioner_Secondary_End_Angle"] = irr_event['positioner_secondary_end_angle'].value[0]
        except:
            df.at[index, "Positioner_Primary_End_Angle"] = float("nan")
            df.at[index, "Positioner_Secondary_End_Angle"] = float("nan")
        df.reset_index(inplace=True, drop=True)
    return(df)
