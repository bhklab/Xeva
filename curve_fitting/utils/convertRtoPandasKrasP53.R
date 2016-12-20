PDXArray = c() # to hold everything
PDXIds = c() # to hold names of PDX models
DataTypes = list("ControlGrowth", "ControlNames", "CaseGrowth", "CaseNames", 
                 "CollectionDays", "CollectionDaysRegrowth", "ContainsRegrowth",
                 "DrugStartDay", "SourceFile")

for (i in 1:length(PharmaData)) {
  PDXIds[i] <- names(PharmaData[i])
}

# create array of all PDX names
for (i in 1:length(PDXIds)) {
  for (j in 1:length(DataTypes)) {
    assign(paste(PDXIds[i], DataTypes[j], sep="_"), c())
  }
}

for (i in 1:length(PDXIds)) {
  for (j in 1:length(DataTypes)) {
    if (DataTypes[j] == "ControlGrowth") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$control_data)
    }
    else if (DataTypes[j] == "CaseGrowth") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$case_data)
    }
    else if (DataTypes[j] == "ControlNames") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$control_names)
    }
    else if (DataTypes[j] == "CaseNames") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$case_names)
    }
    else if (DataTypes[j] == "CollectionDays") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$collection_days)
    }
    else if (DataTypes[j] == "ContainsRegrowth") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$contains_regrowth_phase)
    }
    else if (DataTypes[j] == "DrugStartDay") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$drug_start_day)
    }
    else if (DataTypes[j] == "CaseGrowth") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$case_data)
    }
    else if (DataTypes[j] == "SourceFile") {
      assign(paste(PDXIds[i], DataTypes[j], sep="_"), PharmaData[[i]]$source_file)
    }
  }
}