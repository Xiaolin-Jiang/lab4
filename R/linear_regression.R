#' calculate the statistics using ordinary least squares
#'
#' @param formula .
#' @param data A dataframe.
#' @return the statistics using ordinary least squares.
#' @examples
#' linreg(formula, data)
#' @import methods
#' @importFrom ggplot2 ggplot
#' @importFrom plyr is.formula
#'
#' @export linreg
#' @exportClass linreg

linreg <- setRefClass("linreg",
                      fields = list(beta = "matrix",
                                    y_hat = "matrix",
                                    e_hat = "matrix",
                                    df = "numeric",
                                    sigma_square = "numeric",
                                    variance_beta = "matrix",
                                    t_beta = "matrix",
                                    Form = "formula",
                                    Formulaf = "character",
                                    data = "data.frame",
                                    df_name = "character",
                                    X = "matrix",y = "numeric"),
                      methods = list(
                        initialize = function(Form,data){
                          if(!plyr::is.formula(Form)|!is.data.frame(data))
                            stop()
                          Formulaf <<- Reduce(paste, deparse(Form))
                          df_name <<-  deparse(substitute(data))
                          Form <<-Form
                          data <<- data
                          X <<- model.matrix(Form,data)
                          y_name = all.vars(Form)[1]
                          y <<- data[,y_name]
                          beta <<- solve(t(X) %*% X) %*% t(X) %*% y
                          y_hat <<- X %*% beta
                          e_hat <<- y - y_hat
                          df <<- length(y) - length(all.vars(Form))
                          sigma_square <<- as.numeric(t(e_hat)%*%e_hat)/df
                          variance_beta <<- sigma_square*solve(t(X) %*% X)
                          t_beta <<- beta/sqrt(diag(variance_beta))
                        },

                        show = function() {
                          cat("Call:\n")
                          cat("linreg(formula =",Formulaf,"data=", df_name,")\n\n") #lm(formula = Petal.Length ~ Species, data = iris)
                          cat("Coefficients:\n")
                          #cat(beta)
                          base::print(t(beta))

                        },
                        resid = function(){
                          return(e_hat)
                        },
                        pred = function(){
                          return(y_hat)
                        },

                        coef = function(){
                          coefs = list(beta, y_hat, e_hat, df, sigma_square, variance_beta, t_beta)
                          names(coefs) = c("Regressions coefficients", "fitted values", "residuals", "degrees of freedom",
                                           "residual variance", "variance of the regression coefficients",
                                           "t-values for each coefficient")
                          return(coefs)
                        },

                        summary = function(){
                          # standard error, t-value and p-value as well as the estimate of σˆ and the degrees of freedom in the model.
                          sum_coefs = list(e_hat, df, sigma_square, variance_beta, t_beta)
                          names(sum_coefs) = c("residuals", "degrees of freedom",
                                           "residual variance", "variance of the regression coefficients",
                                           "t-values for each coefficient")
                          base::print(sum_coefs)
                        },

                        plot = function(){
                          df1 = as.data.frame(cbind(y_hat, e_hat))
                          p1 <- ggplot2::ggplot(df1) +
                            ggplot2::aes(x=df1[,1],y=df1[,2])+
                            ggplot2::geom_point(size=2,shape=21)+
                            ggplot2::labs(x=paste0("Fitted values \n linreg(", Formulaf, ")"),
                                          y = "Residuals") +
                            ggplot2::labs(title = "Residuals vs Fitted")+
                            ggplot2::geom_hline(yintercept=0, linetype=3)+
                            ggplot2::stat_summary(fun = median, geom = "line", color = "red")+
                            ggplot2::theme_test()+
                            ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5))


                          df2 = as.data.frame(cbind(y_hat, sqrt(abs(e_hat-mean(e_hat))/sd(e_hat))))
                          p2 <- ggplot2::ggplot(df2) +
                            ggplot2::aes(x=df2[,1],y=df2[,2])+
                            ggplot2::geom_point(size=2,shape=21)+
                            ggplot2::labs(x=paste0("Fitted values \n linreg(", Formulaf, ")"),
                                          y = "sqrt(|Standardized Residuals|))") +
                            ggplot2::labs(title = "Scale−Location")+
                            ggplot2::stat_summary(fun = median, geom = "line", color = "red")+
                            ggplot2::theme_test()+
                            ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5))

                          base::print(p1)
                          base::print(p2)
                        }


                      )
)


