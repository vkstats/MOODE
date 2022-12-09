mood.example <- settings(K = 3, Klev = 3, Nruns = 40)
print.settings(mood.example) #info
mood <- Search(mood.example) #obtain mood for the given settings
mood