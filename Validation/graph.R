library(tidyverse)
library(RColorBrewer)
library(cowplot)

data <- read_csv("Validation.csv")
data <- drop_na(data, percentage_accuracy)
data$percentage_accuracy <- as.numeric(data$percentage_accuracy)

SP1 <- subset(data, Setup == "SP1")
SP2 <- subset(data, Setup == "SP2")
SP3 <- subset(data, Setup == "SP3")
SP <- rbind(SP1, SP2, SP3)
SP_Siemens <- subset(SP, System == "Siemens ARTIS Zee biplane")
SP_GE <- subset(SP, System == "GE Innova IGS 540")
SP_Canon <- subset(SP, System == "Canon Infinix CF-I biplane")
SP_Philips <- subset(SP, System == "Philips Allura Xper")
rm(SP1, SP2, SP3)
colors <- scale_fill_brewer(palette = "Dark2")
axis_ticks = c(0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2)

p1 <- ggplot(SP_Siemens, aes(x = Setup, y = percentage_accuracy, fill = Calculator, color = "black")) +
  geom_col(position = "dodge") + theme_bw() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits = c(0,1.1), breaks=axis_ticks) +
  ggtitle("Simulated procedures (Siemens)") +
  ylab("% of film measurement") +
  xlab(element_blank()) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_color_manual(values = rep("black", 10)) + 
  theme(legend.title=element_blank()) + 
  guides(color = 'none')


p2 <- ggplot(SP_GE, aes(x = Setup, y = percentage_accuracy, fill = Calculator, color = Calculator)) +
  geom_col(position = "dodge") + theme_bw() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits = c(0,1.4), breaks=axis_ticks) +
  ggtitle("Simulated procedures (GE)") +
  ylab("% of film measurement") +
  xlab(element_blank()) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_color_manual(values = rep("black", 10)) + 
  theme(legend.title=element_blank()) + 
  guides(color = 'none')


p3 <- ggplot(SP_Canon, aes(x = Setup, y = percentage_accuracy, fill = Calculator, color = Calculator)) +
  geom_col(position = "dodge") + theme_bw() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits = c(0,1.2), breaks=axis_ticks) +
  ggtitle("Simulated procedures (Canon)") +
  ylab("% of film measurement") +
  xlab(element_blank()) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_color_manual(values = rep("black", 10)) + 
  theme(legend.title=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank()) + 
  guides(color = 'none')

p4 <- ggplot(SP_Philips, aes(x = Setup, y = percentage_accuracy, fill = Calculator, color = Calculator)) +
  geom_col(position = "dodge") + theme_bw() + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0), limits = c(0,1.7), breaks=axis_ticks) +
  ggtitle("Simulated procedures (Philips)") +
  ylab("% of film measurement") +
  xlab(element_blank()) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) + 
  scale_fill_brewer(palette = "Paired") + 
  scale_color_manual(values = rep("black", 10)) + 
  theme(legend.title=element_blank()) + 
  guides(color = 'none')

plot_grid(p1, p2, p3, p4, labels = NULL)