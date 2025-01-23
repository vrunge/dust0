"
library(dust)

y = matrix(rnorm(50), 2)

dust = dust.object.MD(nb_r = 1)
dust$append(y)
dust$update_partition()
dust$get_partition()
dust$get_info()
"
