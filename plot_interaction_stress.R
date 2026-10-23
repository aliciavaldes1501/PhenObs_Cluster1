grid.arrange(
  ggplot() +
    geom_ribbon(data = data.frame(ggpredict(model_threat_mean_slope_all_fl_zi_50_stress,
                                            type = "fixed",
                                            terms=c("mean_slope_abs_fl[all]", "stress_tol[0:96, by=24]"))),
                aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group),
                alpha = 0.25) +
    geom_line(data = data.frame(ggpredict(model_threat_mean_slope_all_fl_zi_50_stress,
                                          type = "fixed",
                                          terms=c("mean_slope_abs_fl[all]", "stress_tol[0:96, by=24]"))),
              aes(x = x, y = predicted, color = group), linewidth = 1) +
    # geom_point(data = sens_per_endanger %>%
    #              filter(slope_var_fl < 50 &
    #                       prop_endangerment > 0),
    #            aes(x = mean_slope_abs_fl, y = prop_endangerment), alpha = 0.6) +
    scale_color_viridis_d(direction = -1) + scale_fill_viridis_d(direction = -1) +
    theme_bw() + ggtitle("p = 0.012") +
    labs(color = "Stress\ntolerance", fill = "Stress\ntolerance") +
    xlab("Sensitivity (absolute value) (flowering)") +
    ylab("Predicted probability of endangerment"),
  ggplot() +
    geom_ribbon(data = data.frame(ggpredict(model_threat_mean_slope_all_fl_zi_50_stress,
                                            type = "zi_prob",
                                            terms=c("mean_slope_abs_fl[all]", "stress_tol[0:96, by=24]"))),
                aes(x = x, y = 1-predicted, ymin = 1-conf.low,
                    ymax = 1-conf.high, fill = group),
                alpha = 0.25) +
    geom_line(data = data.frame(ggpredict(model_threat_mean_slope_all_fl_zi_50_stress,
                                          type = "zi_prob",
                                          terms=c("mean_slope_abs_fl[all]", "stress_tol[0:96, by=24]"))),
              aes(x = x, y = 1-predicted, color = group), linewidth = 1) +
    # geom_point(data = sens_per_endanger %>%
    #              filter(slope_var_fl < 50) %>%
    #              mutate(endangerment_01 = ifelse(prop_endangerment == 0, 0, 1)),
    #            aes(x = mean_slope_abs_fl, y = endangerment_01), alpha = 0.6) +
    scale_color_viridis_d(direction = -1) + scale_fill_viridis_d(direction = -1) +
    theme_bw() + ggtitle("p = 0.127") +
    labs(color = "Stress\ntolerance", fill = "Stress\ntolerance") +
    xlab("Sensitivity (absolute value) (flowering)") +
    ylab("Predicted probability of being endangered"),
  ncol = 2)