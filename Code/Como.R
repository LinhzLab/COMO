
COMO_algo <- function(x,y,wy,Q,beta_inital,rho_inital){
  rho_f <- function(parameter){
    xp = x%*%beta_inital
    expression =(1/(n^2*(n-1)*(n-1)))*(t(y-parameter*wy-xp)%*%Q%*%(y-parameter*wy-xp))^2
  }
  eval_grad_f <- function(parameter) {
    xp = x%*%beta_inital
    y_Q_y = t(y-parameter*wy-xp)%*%Q%*%(y-parameter*wy-xp)
    return((-2/(n^2*(n-1)*(n-1)))*y_Q_y*t(wy)%*%Q%*%(y-parameter*wy-xp))
  }
  solution = optim(rho_inital,fn = rho_f,gr =eval_grad_f,method = 'BFGS')$par
  rho_1 = solution[1]
  beta_f <- function(parameter){
    xp = x%*%parameter
    expression = 1/n*t(y-rho_1*wy-xp)%*%(y-rho_1*wy-xp)+(1/(n^2*(n-1)^2))*(t(y-rho_1*wy-xp)%*%Q%*%(y-rho_1*wy-xp))^2
    return(expression)
  }
  grad_f <- function(parameter){
    xp = x%*%parameter
    y_Q_y = t(y-rho_1*wy-xp)%*%Q%*%(y-rho_1*wy-xp)
    expression = -2*(1/n)*t(x)%*%(y-rho_1*wy-xp)+(-4/(n^2*(n-1)^2))*y_Q_y[1]*t(x)%*%Q%*%(y-rho_1*wy-xp)
    return(expression)
  }
  solution_2 = optim(beta_inital,fn = beta_f,gr = grad_f,method = 'BFGS')$par
  return(list(rho_1 = rho_1, solution_2=solution_2 ))
}