################### Defining stages ####################
stages = list()
stages[["s1"]] = c() # 1 4-7 pcw Embryonic
stages[["s2a"]] = c("8 pcw","9 pcw") # 2A 8-9 pcw Early prenatal
stages[["s2b"]] = c("12 pcw") # 2B 10-12 pcw Early prenatal
stages[["s3a"]] = c("13 pcw") # 3A 13-15 pcw Early mid-prenatal
stages[["s3b"]] = c("16 pcw","17 pcw") # 3B 16-18 pcw Early mid-prenatal
stages[["s4"]] = c("19 pcw","21 pcw","24 pcw") # 4 19-24 pcw Late mid-prenatal
stages[["s5"]] = c("25 pcw","26 pcw","35 pcw","37 pcw") # 5 25-38 pcw Late prenatal
stages[["s6"]] = c("4 mos") # 6 Birth-5 months Early infancy
stages[["s7"]] = c("10 mos") # 7 6-18 months Late infancy
stages[["s8"]] = c("1 yrs","2 yrs","3 yrs","4 yrs") # 8 19 months-5 yrs Early childhood
stages[["s9"]] = c("8 yrs","11 yrs") # 9 6-11 yrs Late childhood
stages[["s10"]] = c("13 yrs","15 yrs","18 yrs","19 yrs") # 10 12-19 yrs Adolescence
stages[["s11"]] = c("21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs") # 11 20-60+ yrs Adulthood
order.stages <- c("s1", "s2a", "s2b", "s3a", "s3b", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11")
