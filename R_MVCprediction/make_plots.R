library(ggpubr)
library(ggside)

#################################################################
# Plot comparison of the models

# Change labels for the facets
session.labs <- as_labeller(c(
  session1 = 'Session 0,8 m/s', 
  session2 = "Session 0,5 m/s"
))

pwc$y.position <- c(0.0535, 0.0682, 0.0654, 0.0712, 0.0668,
                    0.0545, 0.0748, 0.0768, 0.072, 0.0788, 0.0734)

# Plot the results
results.models %>% 
  ggplot(aes(x = model, y = accuracy, color = model)) + 
  geom_boxplot() +     # Add boxplots
  geom_jitter(         # Add the points in front of boxplots
    size = .3,
    alpha = .5
  ) +      
  facet_grid(          # One plot by session
    cols = vars(session),
    labeller = labeller(session = session.labs)
  ) +
  geom_ysidedensity(   # Add distribution at the side of the graphs
    aes(fill = model), 
    color = NA, alpha = 0.4
  ) +
  stat_pvalue_manual(  # Add significance bars
    pwc, tip.length = .01, size = 6
  ) +
  xlab("") + ylab("RMSE") +
  labs(
    caption = expression(~bold("****")~italic(p)~"<"~0.0001~";"~bold("***")~italic(p)~"<"~0.001~";"~bold("**")~italic(p)~"<"~0.01~";"~bold("*")~italic(p)~"<"~0.05)
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    strip.text = element_text(size = 30),
    plot.caption = element_text(size = 20, color = "grey50"),
    ggside.panel.grid = element_blank(),
    ggside.panel.background = element_blank()
  )
  

plot_grid(nrow = 3, aucRMS.VaMe, IMNF.VaLa, aucRMS.VaLa)

plot_grid(nrow = 4, aucRMS.VaMe, aucRMS.ReFe, aucRMS.VaLa, IMNF.ExLo)
#################################################################
# Plot variable importance

# Prepare the trees
tree.session1 <- rpart(formula = estimate.2, data = session.1, 
                       cp = 0.001, minsplit = 20)
rmse.tree.sesion.1 <- 
  sqrt(mean((predict(tree.session1, session.1) - session.1$MVC)^2))

tree.session2 <- rpart(formula = estimate.2, data = session.2, 
                       cp = 0.001, minsplit = 20)
rmse.tree.sesion.2 <- 
  sqrt(mean((predict(tree.session2, session.2) - session.2$MVC)^2))

# Create data frame with variable importance for both sessions
variable.importance.session1 <- as.data.frame(
  tree.session1$variable.importance) %>%
  rownames_to_column(var = "variable")
colnames(variable.importance.session1)[2] <- "value"
variable.importance.session1 <- variable.importance.session1[order(
  variable.importance.session1$value, decreasing = T),]

variable.importance.session1$color <- ifelse(
  startsWith(variable.importance.session1$variable, "IMNF"), "1",
  ifelse(
    startsWith(variable.importance.session1$variable, "aucRMS"), "2", 
    ifelse(
      startsWith(variable.importance.session1$variable, "t2maxRMS"), "3", "4"
      )
    )
  )
  

variable.importance.session2 <- as.data.frame(
  tree.session2$variable.importance) %>%
  rownames_to_column(var = "variable")
colnames(variable.importance.session2)[2] <- "value"
variable.importance.session2 <- variable.importance.session2[order(
  variable.importance.session2$value, decreasing = T),]

variable.importance.session2$color <- ifelse(
  startsWith(variable.importance.session2$variable, "IMNF"), "1",
  ifelse(
    startsWith(variable.importance.session2$variable, "aucRMS"), "2", 
    ifelse(
      startsWith(variable.importance.session2$variable, "t2maxRMS"), "3", "4"
    )
  )
)

# Plot the values
plot.first.variable.importance.tree1 <- variable.importance.session1 %>%
  arrange(desc(value)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(variable, -value, decreasing = TRUE), y = value, fill = factor(color))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  ggtitle(bquote('Arbre de régression session 0,8 '~m.s^{-1}),
          subtitle = paste(
            "RMSE : ", 
            round(rmse.tree.sesion.1, digits = 3))
  )+
  scale_fill_manual(values = c("1" = "dodgerblue", "2" = "lightgreen", "3" = "darkgreen", "4" = "firebrick1")) +
  xlab("") + ylab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 42, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot.first.variable.importance.tree2 <- variable.importance.session2 %>%
  arrange(desc(value)) %>%
  head(5) %>%
  ggplot(aes(x = reorder(variable, -value, decreasing = TRUE), y = value, fill = factor(color))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  ggtitle(bquote('Arbre de régression session 0,5 '~m.s^{-1}),
          subtitle = paste(
            "RMSE : ", 
            round(rmse.tree.sesion.2, digits = 3))
  )+
  scale_fill_manual(values = c("1" = "dodgerblue", "2" = "lightgreen", "3" = "darkgreen", "4" = "firebrick1")) +
  xlab("") + ylab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 42, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


# Fix the number of trees at 300 as mentionned before
kNumberTreeForest <- 300

# Create a random forest by session
forest.session1 <- randomForest(estimate.2, data = session.1, 
                                ntree = kNumberTreeForest)
rmse.forest.sesion.1 <- 
  sqrt(mean((predict(forest.session1, session.1) - session.1$MVC)^2))

forest.session2 <- randomForest(estimate.2, data = session.2, 
                                ntree = kNumberTreeForest)
rmse.forest.sesion.2 <- 
  sqrt(mean((predict(forest.session2, session.2) - session.2$MVC)^2))

# Create a data frame with variable importance for each session
variable.importance.forest.session1 <- data.frame(
  variable = rownames(forest.session1$importance),
  value = forest.session1$importance
)
colnames(variable.importance.forest.session1)[2] <- "value"
variable.importance.forest.session1 <- variable.importance.forest.session1[
  order(variable.importance.forest.session1$value, decreasing = T), ]

variable.importance.forest.session1$color <- ifelse(
  startsWith(variable.importance.forest.session1$variable, "IMNF"), "1",
  ifelse(
    startsWith(variable.importance.forest.session1$variable, "aucRMS"), "2", 
    ifelse(
      startsWith(variable.importance.forest.session1$variable, "t2maxRMS"), "3", "4"
    )
  )
)


variable.importance.forest.session2 <- data.frame(
  variable = rownames(forest.session2$importance), 
  value = forest.session2$importance
)
colnames(variable.importance.forest.session2)[2] <- "value"
variable.importance.forest.session2 <- variable.importance.forest.session2[
  order(variable.importance.forest.session2$value, decreasing = T), ]

variable.importance.forest.session2$color <- ifelse(
  startsWith(variable.importance.forest.session2$variable, "IMNF"), "1",
  ifelse(
    startsWith(variable.importance.forest.session2$variable, "aucRMS"), "2", 
    ifelse(
      startsWith(variable.importance.forest.session2$variable, "t2maxRMS"), "3", "4"
    )
  )
)


# Plot the variable importance
plot.variables.forest.session.1 <- variable.importance.forest.session1 %>%
  head(5) %>%
  ggplot(aes(x = reorder(variable, -value, decreasing = TRUE), y = value, fill = factor(color))) +
  geom_bar(stat = "identity") +
  coord_flip () +
  theme_classic() +
  scale_fill_manual(values = c("1" = "dodgerblue", "2" = "lightgreen", "3" = "darkgreen", "4" = "firebrick1")) +
  xlab("") + ylab("") +
  labs(title = bquote('Forêt alétoire session 0,8 '~m.s^{-1}), 
       subtitle = paste(
         "RMSE : ", 
         round(rmse.forest.sesion.1, digits = 3))
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 42, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot.variables.forest.session.2 <- variable.importance.forest.session2 %>%
  head(5) %>%
  ggplot(aes(x = reorder(variable, -value, decreasing = TRUE), y = value, fill = factor(color))) +
  geom_bar(stat = "identity") +
  coord_flip () +
  theme_classic() +
  labs(title = bquote('Forêt alétoire session 0,5 '~m.s^{-1}),
       subtitle = paste(
         "RMSE : ", 
         round(rmse.forest.sesion.2, digits = 3))
  )+
  xlab("") + ylab("") +
  scale_fill_manual(values = c("1" = "dodgerblue", "2" = "lightgreen", "3" = "darkgreen", "4" = "firebrick1")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 42, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

plot_grid(nrow = 2, ncol = 2, 
          plot.first.variable.importance.tree1, 
          plot.first.variable.importance.tree2,
          plot.variables.forest.session.1, 
          plot.variables.forest.session.2)


variable.importance.forest.session1 %>%
  head(5) %>%
  ggplot(aes(x = reorder(variable, -value, decreasing = TRUE), y = value, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("1" = "dodgerblue", "2" = "lightgreen", "3" = "darkgreen", "4" = "firebrick1")) +
  coord_flip () +
  theme_classic() +
  labs(title = bquote('Forêt alétoire session 0,5 '~m.s^{-1}),
       subtitle = paste(
         "RMSE : ", 
         round(rmse.forest.sesion.2, digits = 3))
  )+
  xlab("") + ylab("")
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 42, hjust = 1, size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

  
  
ggpairs(session.1[,c(5:ncol(session.1))])


session.1 %>%
  filter(subjectName != "01PrAy" | setName != 5 | repetition != 7) %>%
  ggbivariate(outcome = "MVC", explanatory = colnames(session.1)[c(21:25)]) +
  # facet_grid(~subjectName) +
  stat_poly_eq(use_label(c("R2", "P")))


rpart.plot(prune(tree.session1, cp = cp.tree.session1),
           type = 5,     # Displays variable name in interior nodes
           extra = 101,  # Displays number and percentage of observation
           space = 0,    # Deletes the space and make the text bigger
           main = "Session 1"
)
    


x %>%  ggplot(aes(
  x = Muscle, y = Coefficient, shape = signif,
  size = signif, color = model, linewidth = signif)) +
  scale_size_manual(values = c(1.2, 1.2, 1.2, .5)) +
  scale_linewidth_manual(values = c(1.2, 1.2, 1.2, .7)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = .7) +
  geom_pointrange(
    aes(ymin = CI_low, ymax = CI_high),
    position = position_dodge2(reverse = TRUE, width = .5)
    # linewidth = 1.2
  ) +
  scale_shape_manual(values = c(21, 13, 16, 4)) +
  coord_flip() +
  xlab(" ") +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(hjust = 1, size = 10),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15),
    legend.key.size = unit(20, "pt")
  ) +
  facet_wrap(Variable ~ session, ncol = 2, scales = "free")         

x %>%  ggplot(aes(
  x = Muscle, y = Coefficient, shape = factor(signif),
  size = signif, color = session:model, linetype = signif)) +
  scale_size_manual(values = c(1.2, 1.2, 1.2, .5)) +
  scale_linetype_manual(values = c( "solid",  "solid",  "solid", "twodash")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_pointrange(
    aes(ymin = CI_low, ymax = CI_high),
    position = position_dodge2(reverse = TRUE, width = .5)
  ) +
  scale_shape_manual(values = c(21, 13, 16, 4)) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45)
  ) +
  facet_wrap(~Variable, ncol = 1, scales = "free")



cp.tree.session1 = tree.session1$cptable[tree.session1$cptable[,2] == 12, "CP"]

# Find the new cp parameter value
cp.tree.session2 = tree.session2$cptable[tree.session2$cptable[,2] == 12, "CP"]

par(mfrow = c(1,2))
# Plot a tree with the new cp
rpart.plot(prune(tree.session1, cp = cp.tree.session1),
           type = 5,     # Displays variable name in interior nodes
           extra = 1,  # Displays number and percentage of observation
           space = 0,  # Deletes the space and make the text bigger
           main = "Session 0,8 m/s",
           cex = 1.2,
           fallen.leaves = FALSE
)

# Plot a tree with the new cp
rpart.plot(prune(tree.session2, cp = cp.tree.session2),
           type = 5,     # Displays variable name in interior nodes
           extra = 1,  # Displays number and percentage of observation
           space = 0,    # Deletes the space and make the text bigger
           main = "Session 0,5 m/s",
           cex = 1.2,
           fallen.leaves = FALSE
)

prune(tree.session1, cp = cp.tree.session1)$variable.importance %>%
  data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  rename(Importance = ".") %>%
  ggplot(aes(x = fct_reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() +
  xlab(" ")

par(mfrow = c(1, 2))
rpart.plot(prune(tree.session1, cp = cp.tree.session1),
           type = 5,     # Displays variable name in interior nodes
           extra = 101,  # Displays number and percentage of observation
           space = 0,    # Deletes the space and make the text bigger
           main = "Session 1"
)


imp.tree1 <- prune(tree.session1, cp = cp.tree.session1)$variable.importance %>%
  data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  rename(Importance = ".")

imp.tree1 <- imp.tree1[order(-imp.tree1$Importance),] 
imp.tree1$rank <- rank(-imp.tree1$Importance, ties.method = "min")
imp.tree1$session <- "A"

imp.tree2 <- prune(tree.session2, cp = cp.tree.session2)$variable.importance %>%
  data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  rename(Importance = ".")

imp.tree2 <- imp.tree2[order(-imp.tree2$Importance),] 
imp.tree2$rank <- rank(-imp.tree2$Importance, ties.method = "min")
imp.tree2$session <- "B"



library(grid)
library(gridExtra)

g1 <- imp.tree1 %>% 
  head (23) %>%
  ggplot(aes(x = rank, y = Importance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Variable, y = Importance ), hjust = 1.1) +
  coord_flip() +
  scale_x_reverse() +  scale_y_reverse() +
  xlab(" ") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank(),
        plot.margin = unit(c(1,-1,1,0), "mm"))

g2 <- imp.tree2 %>% 
  head (23) %>%
  ggplot(aes(x = rank, y = Importance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Variable, y = Importance + 0.001), hjust = 0) +
  coord_flip() +
  ylim(0, 0.30) +
  scale_x_reverse() +
  xlab(" ") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank(),
        plot.margin = unit(c(1,-1,1,0), "mm"))



gg1 <- ggplot_gtable(ggplot_build(g1))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))
gg2 <- ggplot_gtable(ggplot_build(g2))
grid.arrange(gg1, gg2, ncol = 2)



imp.tree1 <- forest.session1$importance %>%
  data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  rename(Importance = "IncNodePurity")

imp.tree1 <- imp.tree1[order(-imp.tree1$Importance),] 
imp.tree1$rank <- rank(-imp.tree1$Importance, ties.method = "min")
imp.tree1$session <- "A"

imp.tree2 <- forest.session2$importance %>%
  data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  rename(Importance = "IncNodePurity")

imp.tree2 <- imp.tree2[order(-imp.tree2$Importance),] 
imp.tree2$rank <- rank(-imp.tree2$Importance, ties.method = "min")
imp.tree2$session <- "B"



library(grid)
library(gridExtra)

g1 <- imp.tree1 %>% 
  ggplot(aes(x = rank, y = Importance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Variable, y = Importance ), hjust = 1.1) +
  coord_flip() +
  scale_x_reverse() +  scale_y_reverse() +
  ggtitle("Session 0,8 m/s") +
  xlab(" ") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank(),
        plot.margin = unit(c(1,-1,1,0), "mm"),
        plot.title = element_text(hjust = .5))

g2 <- imp.tree2 %>% 
  ggplot(aes(x = rank, y = Importance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Variable, y = Importance + 0.001), hjust = 0) +
  coord_flip() +
  ylim(0, 0.15) +
  ggtitle("Session 0,5 m/s") +
  scale_x_reverse() +
  xlab(" ") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank(),
        plot.title = element_text(hjust = .5),
        plot.margin = unit(c(1,-1,1,0), "mm"))



gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
grid.arrange(gg1, gg2, ncol = 2)



