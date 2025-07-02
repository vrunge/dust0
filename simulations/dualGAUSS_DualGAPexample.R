# Define the objective function (to MINIMIZE, so we negate it)
f <- function(x) {
  x1 <- x[1]
  x2 <- x[2]
  num <- (x1 - x2)^2
  den <- 1 - 2*x1 - x2
  if (den <= 0) return(Inf)  # enforce domain
  -(-1/2 * num / den + 3/2 * (1 + x1 + x2))  # negative for maximization
}

# Initial guess (must satisfy constraints)
x0 <- c(0.3, 0.3)

# Constraint matrix A and vector b for Ax >= b
# Constraints:
#  x1 > 0   →   x1 >= ε  →   1·x1 + 0·x2 ≥ ε
#  x2 > 0   →   x2 >= ε  →   0·x1 + 1·x2 ≥ ε
#  1 - 2x1 - x2 > 0   →   2x1 + x2 ≤ 1  →  -2x1 - x2 ≥ -1
A <- rbind(
  c( 1,  0),  # x1 ≥ ε
  c( 0,  1),  # x2 ≥ ε
  c(-2, -1)   # 2x1 + x2 ≤ 1 → -2x1 - x2 ≥ -1
)

b <- c(0, 0, -1)

# Run optimization
result <- constrOptim(theta = x0, f = f, grad = NULL, ui = A, ci = b)

# Output result
cat("Maximizer:\n")
print(result$par)
cat("Maximum value:\n")
print(-result$value)  # because we minimized the negative

