modelsclimate <- list(
  . ~ . + FFP,
  . ~ . + MAP,
  . ~ . + CMD,
  . ~ . + TD,
  . ~ . + EMT,
  . ~ . + TD + FFP + EMT,
  . ~ . + MAP + CMD,
  . ~ . + TD + FFP + EMT + CMD,
  . ~ . + MAP + FFP + TD + CMD + EMT)
