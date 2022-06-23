preprocess_URs <- function(datasets,activity,Datasets_inf, Datasets_Noninf){
  
  #file names:
  datasets = datasets.drop(21)
  datasets.index = np.array(range(len(datasets)))
  datasets = datasets[0]
  
  
  #if disease is active or not
  activity = activity.drop(21)
  activity.index = np.array(range(len(activity)))
  activity = activity[0]
  
  #Gene translation
  
  #If UR is predicted to be UR in a given dataset for SP1.6
  Datasets_Inf.index = URs
  Datasets_Inf.columns = datasets[activity == 'yes']
  Datasets_Noninf.index = URs
  Datasets_Noninf.columns = datasets[activity == 'no']
  return(list(Datasets_Inf, Datasets_Noninf))
}
