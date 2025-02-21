source("src/utils/get_mse.R")
owrf <- function(x_train, y_train, x_test, y_test, n_tree, norm_param) {
    n_train <- nrow(x_train)
    rf_reg <- randomForest::randomForest(
        x = x_train, y = y_train, mtry = floor(ncol(x_train) / 3), ntree = n_tree,
        replace = TRUE, nodesize = floor(sqrt(n_train)), keep.inbag = TRUE, keep.forest = TRUE
    )

    pred_train <- as.matrix(predict(rf_reg, x_train))
    pred_train_matrix <- predict(rf_reg, x_train, predict.all = TRUE)$individual
    pred_test_matrix <- predict(rf_reg, x_test, predict.all = TRUE)$individual

    num_resample <- rf_reg$inbag
    node_matrix <- matrix(attr(predict(rf_reg, x_train, nodes = TRUE), "nodes"),
        nrow = n_train, ncol = n_tree
    )
    count_in_same_node <- sapply(1:n_tree, function(j) {
        column <- node_matrix[, j]
        unique_elements <- unique(column)

        count_in_same_node_this_col <- sapply(unique_elements, function(x) {
            this_idx <- which(column == x)
            this_count <- sum(num_resample[this_idx, j])
            return(this_count)
        })

        sapply(column, function(x) {
            count_in_same_node_this_col[which(x == unique_elements)]
        })
    })

    m <- num_resample / count_in_same_node

    ## 2 steps OWRF:
    y_hat <- t(pred_train_matrix)
    c0 <- 2 * y_hat %*% t(y_hat) + diag(10^-6, n_tree, n_tree)
    sc0 <- norm(c0, "2")
    c0 <- c0 / sc0
    y_pred_train <- predict(rf_reg, x_train)
    sigma_square_hat_c1 <- ((norm(y_train - pred_train, "2"))**2) / n_train
    d0 <- as.vector((-2 * y_hat %*% y_train + t(m) %*% matrix(sigma_square_hat_c1, nrow = n_train, ncol = 1)) / sc0)
    lb0 <- matrix(0, nrow = n_tree, ncol = 1)
    ub0 <- matrix(1, nrow = n_tree, ncol = 1)
    aeq0 <- matrix(1, nrow = 1, ncol = n_tree)
    beq0 <- matrix(1, nrow = 1, ncol = 1)

    w0_star <- matrix(pracma::quadprog(C = c0, d = d0, Aeq = aeq0, beq = beq0, lb = lb0, ub = ub0)$xmin)
    pred_train_c1_step1 <- pred_train_matrix %*% w0_star
    residual_step1 <- as.matrix(y_train - pred_train_c1_step1)
    ## step2:
    c0_ <- 2 * y_hat %*% t(y_hat) + diag(10^-6, n_tree, n_tree)
    sc0_ <- norm(c0_, "2")
    c0_ <- c0_ / sc0_
    d0_ <- as.vector((-2 * y_hat %*% y_train + t(m) %*% residual_step1**2) / sc0_)
    lb0_ <- matrix(0, nrow = n_tree, ncol = 1)
    ub0_ <- matrix(1, nrow = n_tree, ncol = 1)
    aeq0_ <- matrix(1, nrow = 1, ncol = n_tree)
    beq0_ <- matrix(1, nrow = 1, ncol = 1)

    #  Get weights
    w_star_2steps <- matrix(pracma::quadprog(C = c0_, d = d0_, Aeq = aeq0_, beq = beq0_, lb = lb0_, ub = ub0_)$xmin)

    #  Get predictions
    y_pred <- pred_test_matrix %*% w_star_2steps
    y_pred <- (y_pred * norm_param$sd) + norm_param$mean
    mse_owrf <- mse(y_pred, y_test)
    return(mse_owrf)
}
