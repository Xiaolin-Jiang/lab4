#' calculate the statistics using ordinary least squares
#'
#' @field beta: regression coefficient
#' @field y_hat: fitted value
#' @field e_hat : residuals (e with a hat)
#' @field df : degree of fredom
#' @field sigma_square: variance of residuals
#' @field variance_beta: variance of regression coefficient
#' @field t_beta: t value of regression coefficient
#' @field Form: input formula
#' @field data : input dataframe(the data reference of formula)
#' @field dfname: name of dataframe
#' @field Formulaf : string of formula(without bracket)
#' @field X: created by model.matrix()
#' @field y: true value in dataframe not estimated value
#'
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

                        print = function() {
                          #cat("Call:\n")
                          string1 <- paste0('linreg(formula = ',Formulaf,", data = ", df_name,")") #lm(formula = Petal.Length ~ Species, data = iris)
                          #cat("Coefficients:\n")
                          #cat(beta)
                          base::print(string1)
                          base::print(beta[,1])
                          # should_print <- paste0('linreg(formula = ',Formulaf)
                          # should_print <- paste0(should_print,', data = ')
                          # should_print <- paste0(should_print,df_name)
                          # should_print <- paste0(should_print,')')
                          # base::print(should_print)
                          # base::print(beta[,1])

                        },
                        resid = function(){
                          return(e_hat)
                        },
                        pred = function(){
                          return(y_hat)
                        },

                        coef = function(){
                           return(beta[,1])
                        },

                        summary = function(){
                          p_val = 2*pt(abs(t_beta),df = df,lower.tail = FALSE)
                          std_e <- round((sqrt(diag(variance_beta))),3)
                          #std_e <- sqrt(diag(variance_beta))
                          df3 <- as.data.frame(matrix(c(round(beta,3), std_e, round(t_beta,3), format(p_val, scientific = 2)), ncol = 4))
                          for(i in 1:ncol(X)){
                            df3[i,5]<-"***"
                          }
                          colnames(df3) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)","")
                          rownames(df3) <- rownames(beta)
                          base::print(df3)
                          base::print(paste("Residual standard error:", round(sqrt(sigma_square),3), "on", df, "degrees of freedom"))
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
                            ggplot2::labs(title = "Scale-Location")+
                            ggplot2::stat_summary(fun = median, geom = "line", color = "red")+
                            ggplot2::theme_test()+
                            ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5))

                          base::print(p1)
                          base::print(p2)
                        }


                      )
)


