modelsclimate <- list(
  . ~ . + FFP,
  . ~ . + MAP,
  . ~ . + CMD,
  . ~ . + TD,
  . ~ . + TD + FFP,
  . ~ . + MAP + CMD,
  . ~ . + TD + FFP + CMD,
  . ~ . + MAP + FFP + TD + CMD)
