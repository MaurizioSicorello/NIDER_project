
# create matrix with known correlation
x <- rnorm(100)
y <- rnorm(100)

df <- cbind(x,y)
df_chol <- solve(chol(var(df)))
df_ind <- df%*%df_chol

corrMat <- matrix(c(1,0.15,0.15,1), nrow = 2)  
df_new <- df_ind%*%chol(corrMat)

test <- cor.test(df_new[,1], df_new[,2], "greater")


t_value <- test$statistic

R2 = t_value^2 / (t_value^2 + 100 - 2)
sqrt(R2)


