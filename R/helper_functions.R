
# Assign walleye length categories
assign_wae_psd <- function(data){
  
  ifelse((data>=250)&(data<380), "S-Q",
         ifelse((data>=380)&(data<510), "Q-P",
                ifelse((data>=510)&(data<630), "P-M",
                       ifelse((data>=630)&(data<760), "M-T",
                              ifelse(data>=760, ">T", "SS")))))
}

# Helper for assinging length classes

round_down <- function(x,to=10)
{
  to*(x %/% to)
}