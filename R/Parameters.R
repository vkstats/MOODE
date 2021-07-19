

# Parameters=function(K,index)
# {
#   if(index==1){#quadratic
#     Parameters=c(rep(1,K),0,rep(1,K-1), rep(1,K*(K-1)/2))
# terms: factors, factors^2, 2-fi
#   }
#   if(index==2){#reduced cubic
#     Parameters=c(rep(1,K),0,rep(1,K-1), rep(1,K*(K-1)/2),rep(1,K*(K-1)), rep(1, choose(K,3)))
#     # terms: factors, factors^2, 2-fi, factors^2*factors, 3-fi
#   }
#   return(Parameters)
# }

