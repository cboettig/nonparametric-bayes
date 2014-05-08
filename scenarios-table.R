scenarios <- melt(parameters)

ids <- c(rep(1:3, 3), rep(4:6, 3))

scenarios <- cbind("scenario"= ids, scenarios)
names(scenarios) <- c("scenario", "value", "parameter", "model")
scenarios <- scenarios[c("scenario", "model", "parameter", "value")]

Allen <- scenarios$model == "Allen"
Myers <- scenarios$model == "Myers"

allen_ids <- scenarios[Allen, "parameter"]
myers_ids <- scenarios[Myers, "parameter"]

allen_names <- names(parameters$Allen[[1]])
myers_names <- names(parameters$Myers[[1]])

scenarios[Allen, "parameter"] <- allen_names[allen_ids]
scenarios[Myers, "parameter"] <- myers_names[myers_ids]


scenarios
