def read_RDSR_acc_dose_data(filename):
    import rdsr_navigator as nav
    import pandas as pd
    from pydicom import dcmread
    rdsr_obj = nav.read_file(filename)
    dicom = dcmread(filename, force=True)
    df = pd.DataFrame({
      "Study_Date":pd.Series([],dtype="string"),
      "Study_Description":pd.Series([],dtype="string"),
      "Patient_Name":pd.Series([],dtype="string"),
      "Patient_ID":pd.Series([],dtype="string"),
      "Gender":pd.Series([],dtype="string"),
      "Birth_Date":pd.Series([],dtype="string"),
      "Performing_Physician":pd.Series([],dtype="string"),
      "Referring_Physician":pd.Series([],dtype="string"),
      "Manufacturer":pd.Series([],dtype="string"),
      "Model":pd.Series([],dtype="string"),
      "Serial_Number":pd.Series([],dtype="string"),
      "Fluoro_Dose_RP_Total":pd.Series([],dtype="float64"),
      "Fluoro_Dose_Area_Product_Total":pd.Series([],dtype="float64"),
      "Total_Fluoro_Time":pd.Series([],dtype="float64"),
      "Acquisition_Dose_RP_Total":pd.Series([],dtype="float64"),
      "Acquisition_Dose_Area_Product_Total":pd.Series([],dtype="float64"),
      "Total_Acquisition_Time":pd.Series([],dtype="float64"),
      "Height_Of_System":pd.Series([],dtype="float64"),})
    
    for index, acc_dose in enumerate(rdsr_obj.get_all('accumulated_x-ray_dose_data')):
        df.at[index, "Fluoro_Dose_RP_Total"] = acc_dose['113728:DCM'].value[0]
        df.at[index, "Fluoro_Dose_Area_Product_Total"] = acc_dose['113726:DCM'].value[0]*10000
        df.at[index, "Total_Fluoro_Time"] = acc_dose['113730:DCM'].value[0]
        df.at[index, "Acquisition_Dose_RP_Total"] = acc_dose['113729:DCM'].value[0]
        df.at[index, "Acquisition_Dose_Area_Product_Total"] = acc_dose['113727:DCM'].value[0]*10000
        df.at[index, "Total_Acquisition_Time"] = acc_dose['113855:DCM'].value[0]
        try:
            df.at[index, "Height_Of_System"] = acc_dose['height_of_system'].value[0]
        except:
            df.at[index, "Height_Of_System"] = float("nan")
        df.reset_index(inplace=True, drop=True)
    
    for index in range(0,len(df)):
        df.loc[index, 'Study_Date'] = dicom.StudyDate
        try:
          df.loc[index, "Study_Description"] = dicom.StudyDescription
        except:
          df.loc[index, "Study_Description"] = "Not provided in RDSR"
        df.loc[index, 'Patient_Name'] = str(dicom.PatientName).replace("^", " ")
        df.loc[index, 'Patient_ID'] = dicom.PatientID
        df.loc[index, 'Gender'] = dicom.PatientSex
        df.loc[index, 'Birth_Date'] = dicom.PatientBirthDate
        try:
          df.loc[index, 'Performing_Physician'] = str(dicom.PerformingPhysicianName).replace("^", " ")
        except:
          df.loc[index, 'Performing_Physician'] = "Not provided in RDSR"
        try:
          df.loc[index, 'Referring_Physician'] = dicom.ReferringPhysicianName.replace("^", " ")
        except:
          df.loc[index, 'Referring_Physician'] = "Not provided in RDSR"
        try:
          df.loc[index, 'Manufacturer'] = dicom.Manufacturer
        except:
          df.loc[index, 'Manufacturer'] = "Not provided in RDSR"
        try:
          df.loc[index, 'Model'] = dicom.ManufacturerModelName
        except:
          df.loc[index, 'Model'] = "Not provided in RDSR"
        try:
          df.loc[index, 'Serial_Number'] = dicom.DeviceSerialNumber
        except:
          df.loc[index, 'Serial_Number'] = "Not provided in RDSR"
    return(df)
