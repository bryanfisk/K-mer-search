args = commandArgs(TRUE)

x = as.double(args[1])
t = as.double(args[2])
r = as.double(args[3])

poisson.test(x, t, r)