set.seed(123)
if (!requireNamespace("e1071", quietly = TRUE)) install.packages("e1071")
library(e1071)

# 数据预处理
data(iris)
iris_sub <- iris[iris$Species %in% c("versicolor", "virginica"), ]
y_raw <- ifelse(iris_sub$Species == "versicolor", 1, -1)
X_raw <- scale(as.matrix(iris_sub[, 1:4]), center = TRUE, scale = TRUE)

# 8:2比例划分训练集与测试集
n_total <- nrow(X_raw)
train_idx <- sample(seq_len(n_total), size = as.integer(0.8 * n_total))

X_tr <- X_raw[train_idx, ]
y_tr <- y_raw[train_idx]
X_te <- X_raw[-train_idx, ]
y_te <- y_raw[-train_idx]

# 对偶二次型 Q 矩阵
Q_tr <- (y_tr %*% t(y_tr)) * (X_tr %*% t(X_tr))

# 双重约束投影算子
project_svm_constraints <- function(v, y, C) {
  v_box <- pmin(pmax(v, 0), C)
  f_theta <- function(theta) {
    sum(pmin(pmax(v - theta * y, 0), C) * y)
  }
  search_bound <- max(100 * C, 1000)
  root_theta <- tryCatch({
    if (f_theta(-search_bound) * f_theta(search_bound) < 0) {
      uniroot(f_theta, c(-search_bound, search_bound), tol = 1e-9)$root
    } else { 0 }
  }, error = function(e) { return(0) })
  return(pmin(pmax(v - root_theta * y, 0), C))
}

# 算法1-投影梯度法
pgd_svm <- function(Q, y, C, tol = 1e-6, maxit = 20000) {
  start_time <- Sys.time()
  N <- nrow(Q)
  a_vec <- rep(0, N)
  
  fval_seq <- numeric(maxit + 1)
  fval_seq[1] <- 0.5 * as.numeric(t(a_vec) %*% Q %*% a_vec) - sum(a_vec)
  
  for (iter in 1:maxit) {
    Qa <- Q %*% a_vec
    grad <- as.vector(Qa - 1)
    f_cur <- fval_seq[iter]
    
    step <- 1
    sigma <- 1e-4
    gamma_bt <- 0.5
    for (bt in 1:30) {
      a_try <- project_svm_constraints(a_vec - step * grad, y, C)
      f_try <- 0.5 * sum(a_try * (Q %*% a_try)) - sum(a_try)
      if (f_try <= f_cur + sigma * sum(grad * (a_try - a_vec))) break
      step <- step * gamma_bt
    }
    
    a_new <- project_svm_constraints(a_vec - step * grad, y, C)
    fval_seq[iter + 1] <- 0.5 * sum(a_new * (Q %*% a_new)) - sum(a_new)
    
    gmap_norm <- sqrt(sum((a_vec - a_new)^2)) / step
    if (gmap_norm < tol) { a_vec <- a_new; break }
    a_vec <- a_new
  }
  return(list(alpha = a_vec, iter = iter, fval = fval_seq[1:(iter + 1)], time = as.numeric(Sys.time() - start_time)))
}

# 算法2-近端梯度法ISTA
ista_svm <- function(Q, y, C, tol = 1e-6, maxit = 20000) {
  start_time <- Sys.time()
  N <- nrow(Q)
  a_vec <- rep(0, N)
  
  L <- max(eigen(Q, symmetric = TRUE, only.values = TRUE)$values)
  gamma <- 1 / L
  
  fval_seq <- numeric(maxit + 1)
  fval_seq[1] <- 0.5 * as.numeric(t(a_vec) %*% Q %*% a_vec) - sum(a_vec)
  
  for (iter in 1:maxit) {
    grad <- as.vector(Q %*% a_vec - 1)
    a_old <- a_vec
    
    y_k <- a_old - gamma * grad
    a_new <- project_svm_constraints(y_k, y, C)
    fval_seq[iter + 1] <- 0.5 * sum(a_new * (Q %*% a_new)) - sum(a_new)
    
    gmap <- (a_old - a_new) / gamma
    if (sqrt(sum(gmap^2)) < tol) { a_vec <- a_new; break }
    a_vec <- a_new
  }
  return(list(alpha = a_vec, iter = iter, fval = fval_seq[1:(iter + 1)], time = as.numeric(Sys.time() - start_time)))
}

# 算法3-加速近端梯度法FISTA
fista_svm <- function(Q, y, C, tol = 1e-6, maxit = 20000) {
  start_time <- Sys.time()
  N <- nrow(Q)
  
  L <- max(eigen(Q, symmetric = TRUE, only.values = TRUE)$values)
  gamma <- 1 / L
  
  beta <- rep(0, N) 
  z <- beta          
  t_val <- 1
  
  fval_seq <- numeric(maxit + 1)
  fval_seq[1] <- 0.5 * as.numeric(t(beta) %*% Q %*% beta) - sum(beta)
  
  for (iter in 1:maxit) {
    grad_z <- as.vector(Q %*% z - 1)
    y_k <- z - gamma * grad_z
    
    beta_new <- project_svm_constraints(y_k, y, C)
    fval_seq[iter + 1] <- 0.5 * sum(beta_new * (Q %*% beta_new)) - sum(beta_new)
    
    gmap <- (z - beta_new) / gamma
    if (sqrt(sum(gmap^2)) < tol) { beta <- beta_new; break }
    
    t_next <- (1 + sqrt(1 + 4 * t_val^2)) / 2
    z <- beta_new + ((t_val - 1) / t_next) * (beta_new - beta)
    
    beta <- beta_new
    t_val <- t_next
  }
  return(list(alpha = beta, iter = iter, fval = fval_seq[1:(iter + 1)], time = as.numeric(Sys.time() - start_time)))
}

# 算法3(改进)-带自适应重启机制的改进加速近端梯度法
fista_svm_restart <- function(Q, y, C, tol = 1e-6, maxit = 20000) {
  start_time <- Sys.time()
  N <- nrow(Q)
  
  L <- max(eigen(Q, symmetric = TRUE, only.values = TRUE)$values)
  gamma <- 1 / L
  
  beta <- rep(0, N) 
  z <- beta          
  t_val <- 1
  
  fval_seq <- numeric(maxit + 1)
  fval_seq[1] <- 0.5 * as.numeric(t(beta) %*% Q %*% beta) - sum(beta)
  
  for (iter in 1:maxit) {
    grad_z <- as.vector(Q %*% z - 1)
    y_k <- z - gamma * grad_z
    
    beta_new <- project_svm_constraints(y_k, y, C)
    fval_seq[iter + 1] <- 0.5 * sum(beta_new * (Q %*% beta_new)) - sum(beta_new)
    
    gmap <- (z - beta_new) / gamma
    if (sqrt(sum(gmap^2)) < tol) { beta <- beta_new; break }
    
    # 检测到非单调过冲时清零动量
    if (fval_seq[iter + 1] > fval_seq[iter]) {
      t_val <- 1
      z <- beta_new
    }
    
    t_next <- (1 + sqrt(1 + 4 * t_val^2)) / 2
    z <- beta_new + ((t_val - 1) / t_next) * (beta_new - beta)
    
    beta <- beta_new
    t_val <- t_next
  }
  return(list(alpha = beta, iter = iter, fval = fval_seq[1:(iter + 1)], time = as.numeric(Sys.time() - start_time)))
}

# 辅助函数-自适应参数还原与KKT计算
evaluate_svm_model <- function(res, X_tr, y_tr, X_te, y_te, C, tol = 1e-6) {
  w <- colSums(res$alpha * y_tr * X_tr)
  eps <- tol * 10
  sv_idx <- which(res$alpha > eps & res$alpha < (C - eps))
  if (length(sv_idx) > 0) {
    b <- mean(y_tr[sv_idx] - X_tr[sv_idx, ] %*% w)
  } else {
    b <- mean(y_tr[res$alpha > eps] - X_tr[res$alpha > eps, ] %*% w)
  }
  pred_te <- X_te %*% w + b
  acc <- mean(ifelse(pred_te >= 0, 1, -1) == y_te)
  
  margin_term <- as.vector(y_tr * (X_tr %*% w + b) - 1)
  eq_residual <- sum(res$alpha * y_tr)
  free_kkt <- if(length(sv_idx) > 0) max(abs(res$alpha[sv_idx] * margin_term[sv_idx])) else 0
  bound_kkt <- max(pmax(margin_term[res$alpha >= (C - eps)], 0))
  
  return(list(w = w, b = b, acc = acc, eq = eq_residual, free_kkt = free_kkt, bound_kkt = bound_kkt))
}

# 实证函数
run_comprehensive_empirical <- function() {
  C_default <- 1.0
  tol_default <- 1e-6
  
  # 运行所有手写算法
  res_pgd      <- pgd_svm(Q_tr, y_tr, C = C_default, tol = tol_default)
  res_ista     <- ista_svm(Q_tr, y_tr, C = C_default, tol = tol_default)
  res_fista    <- fista_svm(Q_tr, y_tr, C = C_default, tol = tol_default)
  res_f_rest   <- fista_svm_restart(Q_tr, y_tr, C = C_default, tol = tol_default)
  
  eval_pgd     <- evaluate_svm_model(res_pgd, X_tr, y_tr, X_te, y_te, C_default, tol_default)
  eval_ista    <- evaluate_svm_model(res_ista, X_tr, y_tr, X_te, y_te, C_default, tol_default)
  eval_fista   <- evaluate_svm_model(res_fista, X_tr, y_tr, X_te, y_te, C_default, tol_default)
  eval_f_rest  <- evaluate_svm_model(res_f_rest, X_tr, y_tr, X_te, y_te, C_default, tol_default)
  
  # e1071标准结果
  lib_model <- svm(X_tr, as.factor(y_tr), kernel = "linear", cost = C_default, scale = FALSE)
  w_lib <- as.vector(t(lib_model$coefs) %*% lib_model$SV)
  b_lib <- -lib_model$rho
  pred_lib <- X_te %*% w_lib + b_lib
  acc_lib <- mean(ifelse(pred_lib >= 0, 1, -1) == y_te)
  
# 收敛曲线图
  old_par <- par(bg = "white", fg = "black")
  dev.hold()
  max_len <- max(length(res_ista$fval), length(res_pgd$fval), length(res_fista$fval), length(res_f_rest$fval))
  plot(1:max_len, rep(NA, max_len), xlim = c(1, max_len), ylim = range(res_ista$fval),
       log = "x", xlab = "迭代序号 (k + 1, 对数坐标轴 Log Scale)", ylab = "对偶目标函数值 f(alpha)",
       main = "凸优化算法对偶收敛轨迹深度对比图", type = "n")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  
  lines(1:length(res_ista$fval), res_ista$fval, col = "red", lwd = 1.5)
  lines(1:length(res_pgd$fval), res_pgd$fval, col = "darkgreen", lwd = 1.5)
  lines(1:length(res_fista$fval), res_fista$fval, col = "blue", lwd = 1.5)
  lines(1:length(res_f_rest$fval), res_f_rest$fval, col = "purple", lwd = 2.5, lty = 2)
  
  leg_text <- c(paste("ISTA (", res_ista$iter, "次)"), 
                paste("PGD (", res_pgd$iter, "次)"), 
                paste("经典 FISTA (", res_fista$iter, "次)"),
                paste("改进 FISTA_Restart (", res_f_rest$iter, "次)"))
  max_text_w <- max(strwidth(leg_text, cex = 0.85)) * 1.25
  legend("topright", legend = leg_text, 
         col = c("red", "darkgreen", "blue", "purple"), 
         lwd = c(1.5, 1.5, 1.5, 2.5), 
         lty = c(1, 1, 1, 2), 
         bg = "white",
         cex = 0.85,
         text.width = max_text_w)
  dev.flush()
  par(old_par)

# 决策边界图
  w_2d <- eval_f_rest$w[1:2]
  b_2d <- eval_f_rest$b
  old_par_boundary <- par(bg = "white", fg = "black")
  
  plot(X_tr[, 1], X_tr[, 2], type = "n",
       xlab = "标准化 Sepal.Length (x1)", ylab = "标准化 Sepal.Width (x2)",
       main = "手写凸优化 SVM 决策边界与支持向量几何解析图")
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white", border = NA)
  grid(col = "lightgray", lty = "dotted")
  points(X_tr[y_tr == 1, 1], X_tr[y_tr == 1, 2], col = "blue", pch = 15, cex = 1.1)
  points(X_tr[y_tr == -1, 1], X_tr[y_tr == -1, 2], col = "red", pch = 16, cex = 1.1)
  slope <- -w_2d[1] / w_2d[2]
  intercept_decision <- -b_2d / w_2d[2]
  
  abline(a = intercept_decision, b = slope, col = "black", lwd = 3)
  abline(a = (1 - b_2d) / w_2d[2], b = slope, col = "darkgray", lty = 2, lwd = 1.5)
  abline(a = (-1 - b_2d) / w_2d[2], b = slope, col = "darkgray", lty = 2, lwd = 1.5)
  
  sv_points_idx <- which(res_f_rest$alpha > 1e-5)
  points(X_tr[sv_points_idx, 1], X_tr[sv_points_idx, 2], 
         col = "chartreuse4", lwd = 2, cex = 1.8, pch = 1) 
  # 图例
  legend("bottomleft", 
         legend = c("Versicolor (+1)", "Virginica (-1)", "决策边界", "间隔带边界", "识别的支持向量"),
         col = c("blue", "red", "black", "darkgray", "chartreuse4"),
         pch = c(15, 16, NA, NA, 1),
         lty = c(NA, NA, 1, 2, NA),
         lwd = c(NA, NA, 3, 1.5, 2),
         bg = "white", cex = 0.85)
  par(old_par_boundary)
    
  # 表格1-计算效率对比
  base_metrics <- data.frame(
    Indicator = c("Sepal.Length (w1)", "Sepal.Width (w2)", "Petal.Length (w3)", "Petal.Width (w4)", 
                  "Intercept (b)", "迭代总轮数 (次)", "算法运行耗时 (秒)", "最终对偶目标函数值", 
                  "测试集泛化准确率"),
    Lib_e1071 = c(w_lib[1], w_lib[2], w_lib[3], w_lib[4], b_lib, NA, NA, NA, acc_lib),
    algorithm_PGD = c(eval_pgd$w[1], eval_pgd$w[2], eval_pgd$w[3], eval_pgd$w[4], eval_pgd$b, res_pgd$iter, res_pgd$time, tail(res_pgd$fval,1), eval_pgd$acc),
    algorithm_ISTA = c(eval_ista$w[1], eval_ista$w[2], eval_ista$w[3], eval_ista$w[4], eval_ista$b, res_ista$iter, res_ista$time, tail(res_ista$fval,1), eval_ista$acc),
    algorithm_FISTA = c(eval_fista$w[1], eval_fista$w[2], eval_fista$w[3], eval_fista$w[4], eval_fista$b, res_fista$iter, res_fista$time, tail(res_fista$fval,1), eval_fista$acc),
    algorithm_FISTA_Restart = c(eval_f_rest$w[1], eval_f_rest$w[2], eval_f_rest$w[3], eval_f_rest$w[4], eval_f_rest$b, res_f_rest$iter, res_f_rest$time, tail(res_f_rest$fval,1), eval_f_rest$acc)
  )
  
  # 表格2-收敛精度与分层KKT条件对比
  kkt_metrics <- data.frame(
    Indicator = c("对偶等式约束残差 (sum(alpha*y))", "自由支持向量 KKT 最大残差", 
                  "边界支持向量 KKT 违背量", "w1 绝对误差 (vs e1071, 下同)", "w2 绝对误差", 
                  "w3 绝对误差", "w4 绝对误差", "b 绝对误差"),
    Lib_e1071 = rep(0, 8),
    algorithm_PGD = c(eval_pgd$eq, eval_pgd$free_kkt, eval_pgd$bound_kkt, abs(eval_pgd$w - w_lib), abs(eval_pgd$b - b_lib)),
    algorithm_ISTA = c(eval_ista$eq, eval_ista$free_kkt, eval_ista$bound_kkt, abs(eval_ista$w - w_lib), abs(eval_ista$b - b_lib)),
    algorithm_FISTA = c(eval_fista$eq, eval_fista$free_kkt, eval_fista$bound_kkt, abs(eval_fista$w - w_lib), abs(eval_fista$b - b_lib)),
    algorithm_FISTA_Restart = c(eval_f_rest$eq, eval_f_rest$free_kkt, eval_f_rest$bound_kkt, abs(eval_f_rest$w - w_lib), abs(eval_f_rest$b - b_lib))
  )
  
  # 表格3-参数敏感性分析
  C_test <- 0.1
  res_pgd_c    <- pgd_svm(Q_tr, y_tr, C = C_test, tol = tol_default)
  res_ista_c   <- ista_svm(Q_tr, y_tr, C = C_test, tol = tol_default)
  res_fista_c  <- fista_svm(Q_tr, y_tr, C = C_test, tol = tol_default)
  res_f_rest_c <- fista_svm_restart(Q_tr, y_tr, C = C_test, tol = tol_default)
  
  eval_pgd_c   <- evaluate_svm_model(res_pgd_c, X_tr, y_tr, X_te, y_te, C_test, tol_default)
  eval_ista_c  <- evaluate_svm_model(res_ista_c, X_tr, y_tr, X_te, y_te, C_test, tol_default)
  eval_fista_c <- evaluate_svm_model(res_fista_c, X_tr, y_tr, X_te, y_te, C_test, tol_default)
  eval_f_rest_c<- evaluate_svm_model(res_f_rest_c, X_tr, y_tr, X_te, y_te, C_test, tol_default)
  
  lib_model_c  <- svm(X_tr, as.factor(y_tr), kernel = "linear", cost = C_test, scale = FALSE)
  w_lib_c <- as.vector(t(lib_model_c$coefs) %*% lib_model_c$SV)
  b_lib_c <- -lib_model_c$rho
  
  sensitivity_metrics <- data.frame(
    Indicator = c("C=0.1 下 w1 参数值", "C=0.1 下 w2 参数值", "C=0.1 下 w3 参数值", "C=0.1 下 w4 参数值", 
                  "C=0.1 下截距 b 值", "该惩罚项下收敛迭代轮数 (次)"),
    Lib_e1071 = c(w_lib_c[1], w_lib_c[2], w_lib_c[3], w_lib_c[4], b_lib_c, NA),
    algorithm_PGD = c(eval_pgd_c$w[1], eval_pgd_c$w[2], eval_pgd_c$w[3], eval_pgd_c$w[4], eval_pgd_c$b, res_pgd_c$iter),
    algorithm_ISTA = c(eval_ista_c$w[1], eval_ista_c$w[2], eval_ista_c$w[3], eval_ista_c$w[4], eval_ista_c$b, res_ista_c$iter),
    algorithm_FISTA = c(eval_fista_c$w[1], eval_fista_c$w[2], eval_fista_c$w[3], eval_fista_c$w[4], eval_fista_c$b, res_fista_c$iter),
    algorithm_FISTA_Restart = c(eval_f_rest_c$w[1], eval_f_rest_c$w[2], eval_f_rest_c$w[3], eval_f_rest_c$w[4], eval_f_rest_c$b, res_f_rest_c$iter)
  )
  
  # 表格4-算法稳定性分析
  tol_strict <- 1e-8
  res_pgd_s    <- pgd_svm(Q_tr, y_tr, C = C_default, tol = tol_strict)
  res_ista_s   <- ista_svm(Q_tr, y_tr, C = C_default, tol = tol_strict)
  res_fista_s  <- fista_svm(Q_tr, y_tr, C = C_default, tol = tol_strict)
  res_f_rest_s <- fista_svm_restart(Q_tr, y_tr, C = C_default, tol = tol_strict)
  
  eval_pgd_s    <- evaluate_svm_model(res_pgd_s, X_tr, y_tr, X_te, y_te, C_default, tol_strict)
  eval_ista_s   <- evaluate_svm_model(res_ista_s, X_tr, y_tr, X_te, y_te, C_default, tol_strict)
  eval_fista_s  <- evaluate_svm_model(res_fista_s, X_tr, y_tr, X_te, y_te, C_default, tol_strict)
  eval_f_rest_s <- evaluate_svm_model(res_f_rest_s, X_tr, y_tr, X_te, y_te, C_default, tol_strict)
  
  stability_metrics <- data.frame(
    Indicator = c("高精度1e-8下迭代轮数", "高精度下自由支持向量 KKT 残差", "高精度下 w1 估计值", "高精度下截距 b 估计值"),
    Lib_e1071 = c(NA, 0, w_lib[1], b_lib),
    algorithm_PGD = c(res_pgd_s$iter, eval_pgd_s$free_kkt, eval_pgd_s$w[1], eval_pgd_s$b),
    algorithm_ISTA = c(res_ista_s$iter, eval_ista_s$free_kkt, eval_ista_s$w[1], eval_ista_s$b),
    algorithm_FISTA = c(res_fista_s$iter, eval_fista_s$free_kkt, eval_fista_s$w[1], eval_fista_s$b),
    algorithm_FISTA_Restart = c(res_f_rest_s$iter, eval_f_rest_s$free_kkt, eval_f_rest_s$w[1], eval_f_rest_s$b)
  )
  
  # 输出结论
  cat("一、 收敛曲线图见弹出窗口\n")
  
  cat("二、 计算效率对比\n")
  print(base_metrics, digits = 5, row.names = FALSE)
  cat("\n")
  
  cat("三、 收敛精度与分层 KKT 条件实证\n")
  print(kkt_metrics, digits = 5, row.names = FALSE)
  cat("\n")
  
  cat("四、 参数敏感性分析\n")
  print(sensitivity_metrics, digits = 5, row.names = FALSE)
  cat("\n")
  
  cat("五、 算法稳定性分析\n")
  print(stability_metrics, digits = 5, row.names = FALSE)
  cat("\n")
}

run_comprehensive_empirical()
