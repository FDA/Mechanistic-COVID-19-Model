library(ggplot2)

risknames <- c("low","intermediate","high")
risks <- risks[names(lb)]
df <- tibble(ID = names(lb),lower = lb,avg = med,upper = ub,riskcategories = risknames[risks + 1])

jindex <- order(df$riskcategories,decreasing = TRUE)
df <- df[jindex,]

#jindex2 <- order(df$upper,decreasing=TRUE)
#df <- df[jindex,]

#df$riskcategories <- ordered(x = df$riskcategories, levels = risknames)
df$ID <- ordered(df$ID,levels=df$ID)

myplot <- ggplot(data = df, aes(x = ID, y = avg, ymin = lower, ymax = upper, colour = riskcategories)) +
  geom_point(position = position_dodge(width = 0.2), size = 3.5) +
  geom_errorbar(position = position_dodge(width = 0.2), width = 0.8, size = 1.0) +
  coord_flip() +
  scale_colour_manual(values = c("red", "blue","green")) +
   ylab("Safety Margin") + xlab("Drug Name") + ylim(0,1000)

print(myplot)