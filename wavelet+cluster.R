library(WaveletComp)
library(tidyr)
library(readxl)
library(dtwclust) #для кластеризации, в том числе и временных рядов
library(dplyr)
library(purrr) #для функции map
library(lubridate)
library(ggplot2)
library(tseries) # Проверка временных рядов на стационарность
library(ggfortify)

################################################################################
##########                 Import data                    ######################
setwd("A:/numbers_of_mammals/wavelet+cluster+autocorrelation")

# XX век
XX_both_bank <- as.data.frame(read_excel("data/XX/All_XX.xlsx"))
LB_XX <- as.data.frame(read_excel("data/XX/LB_XX.xlsx"))
RB_XX <- as.data.frame(read_excel("data/XX/RB_XX.xlsx"))

# XXI век
XXI_both_bank <- as.data.frame(read_excel("data/XXI/XXI_all.xlsx"))
LB_XXI <- as.data.frame(read_excel("data/XXI/XXI_LB.xlsx"))
RB_XXI <- as.data.frame(read_excel("data/XXI/XXI_RB.xlsx"))

#########################################################################
########## Объединенный датафрейм для вейвлет-преобразования (между веками нули) 
############
XX_both_bank_long <- as.data.frame(read_excel("data/XX/All_XX.xlsx"))%>%
  pivot_longer(cols = -date, names_to="Spec", values_to = "N")
LB_XX_long <- as.data.frame(read_excel("data/XX/LB_XX.xlsx"))%>%
  pivot_longer(cols = -date, names_to="Spec", values_to = "N")
RB_XX_long <- as.data.frame(read_excel("data/XX/RB_XX.xlsx"))%>%
  pivot_longer(cols = -date, names_to="Spec", values_to = "N")
XXI_both_bank_long <- as.data.frame(read_excel("data/XXI/XXI_all.xlsx"))%>%
  pivot_longer(cols = -date, names_to="Spec", values_to = "N")
LB_XXI_long <- as.data.frame(read_excel("data/XXI/XXI_LB.xlsx"))%>%
  pivot_longer(cols = -date, names_to="Spec", values_to = "N")
RB_XXI_long <- as.data.frame(read_excel("data/XXI/XXI_RB.xlsx"))%>%
  pivot_longer(cols = -date, names_to="Spec", values_to = "N")
zero_between_century <- data.frame(date=c(1995:2007), 
                                   Sor_ara = 0, Sor_cae = 0, Sor_min = 0,
                                   Sor_iso = 0, Sor_rob = 0, Sor_tun = 0, 
                                   Sor_mss = 0, Sor_dap = 0, Neo_fod = 0,
                                   Tal_alt = 0, Cle_rut = 0, Cle_gla = 0,
                                   Cle_rfc = 0, Mic_agr = 0, Mic_oec = 0,
                                   Myo_sch = 0, Mcr_min = 0, Sic_bet = 0,
                                   Community = 0)%>%
  pivot_longer(cols = -date, names_to="Spec", values_to = "N")

both_bank0 <- rbind(XX_both_bank_long, zero_between_century, XXI_both_bank_long)%>%
  pivot_wider(names_from = "Spec", values_from = "N")%>%
  select(-Arv_ter,-Tal_alt, -Mcr_min)
LB0 <- rbind(LB_XX_long, zero_between_century, LB_XXI_long)%>%
  pivot_wider(names_from = "Spec", values_from = "N")%>%
  select(-Arv_ter, -Tal_alt, -Mcr_min)

RB0 <- rbind(RB_XX_long, zero_between_century, RB_XXI_long)%>%
  pivot_wider(names_from = "Spec", values_from = "N")%>%
  select(-Arv_ter, -Tal_alt, -Mcr_min)

rm(LB_XX_long, LB_XXI_long, RB_XX_long, RB_XXI_long, XX_both_bank_long,
   XXI_both_bank_long, zero_between_century)




################################################################################
##########                 Function                         ####################
#подумать над проверками исходного датафрейма
# подумать, почему не сохраняет вывод в переменную

wavelet_cluster <- function(df,rename){
  if (df[1,1] > 1999){
    cat <- "XXI"
  } 
  else {
    cat <- "XX"
  }
  df <- df%>%
    rename_with(~paste0(rename, .x), .cols=-1)%>%
    mutate(date = as.Date(as.character(date), format = "%Y"))
  
  graph_wavelet <- function(df,col, rename){
    wt.image(analyze.wavelet(df, col, upperPeriod = 5), color.key="interval", 
             show.date = T,main = col, legend.params=list(lab="wavelet power levels"))
    
    
    dev.print(png, filename = paste(cat, '/', rename, '/', col, '.png', sep=''), width = 1000, height = 500)
  } #функция для графика вейвлет-спектра
  pwr_wavelet <- function(df,col){
    return(analyze.wavelet(df,col,upperPeriod = 5)$Power)
  } # функция для вейвлет преобразования, на выходе power
  my_cluster <- function(x){
    result_of_cluster <- tsclust(
      x,
      k = 6,                 # запрашиваемое число кластеров
      type = "hierarchical", # тип кластеризации
      distance = "dtw",      # мера расстояния - динамическая трансформации временной шкалы (Dynamic Time Warping)
      seed = 42,
      control = 
        hierarchical_control(method = "ward.D2"), # метод агломерации
      args = 
        tsclust_args(dist = list(window.size = 7)) # размер окна Сакэ-Чиба
    )
    return(result_of_cluster)
  }#запустить кластеризацию
  plot_of_cluster <- function(x, names){
    return(plot(x, xlab = "", sub = "", main = names))
  } #построить дендрограмму
  
  
  
  graph <- map(names(df)[2:ncol(df)], graph_wavelet, df = df, rename = rename) #все графики спектра вейвлет
  names(graph) <- names(df)[2:ncol(df)] # задал нормальные названия
  
  power_wavelet <- map(names(df)[2:ncol(df)], pwr_wavelet, df = df) #power спектры вейвлета для всех переменных
  names(power_wavelet) <- names(df)[2:ncol(df)]
  cluster <- my_cluster(power_wavelet) #результат кластеризации
  dendrogram <- plot_of_cluster(cluster, "Wavelet power") #дендрограмма по спектрам вейвлета
  
  
  numbers <- as.list(df[2:ncol(df)])
  cluster_of_number <- my_cluster(numbers) 
  
  dendrogram_of_numbers <- plot_of_cluster(cluster_of_number, "Кластеризация численности")
  
  return(dendrogram_of_numbers)
  return(cluster)
  return(graph)
  return(dendrogram)
}   # rename - просто добавляет к столбцам индекс берега (какой задашь)
test_autocorrelation <- function(y){
  names <- colnames(y)
  new <- rbind(names,y)
  autocorr <- function(x){ acf(as.numeric(x[2:length(x)]), lag.max = 6, main=x[1])}
  result <- apply(new, 2, autocorr)
  return(result)
} #функци¤ строит функции автокорреляции дл¤ каждого столбца
#при этом называет график автокорреляции именем столбца
# выводит результаты в переменную result
#нужно сделать таблицу с коэффициентами автокоррелляции, уровнем значимости
dominant_structure_for_graph <- function(df){
  community <- df%>%
    select(date,Community)
  df2 <- df%>%
    mutate_at(vars(c(2:ncol(df))), funs(.*100/df[,ncol(df)]))%>%
    mutate_if(is.numeric, round, digits=1)%>%
    pivot_longer(cols = -date, names_to="Spec", values_to = "Dominant")
  df3 <- df%>%
    mutate_if(is.numeric, round, digits=1)%>%
    pivot_longer(cols = -date, names_to="Spec", values_to = "Numbers")
  
  df4 <- full_join(df2,df3, by=c("date", "Spec"))
  df5 <- full_join(df4,community, by = "date")
  return(df5)
} #
# Объединенная таблица: численность и структура доминирования (доля вида в конкретный год),
# при этом входной
# датафрейм должен в первом столбце иметь год, в последнем - общую численность сообщества
plot_numbers_and_dominant <- function(x, n_cols){
  return(
    ggplot(x)+
      geom_col(aes(date, Dominant), fill = "646464")+
      geom_point(aes(date, Numbers))+
      geom_line(aes(date, Numbers))+
      scale_y_continuous(
        "Численность (экз-100 ц/с)",
        sec.axis=sec_axis(~.,name="Доля вида в сообществе, %"))+
      facet_wrap(~Spec, ncol=n_cols)+
      theme_test()
  )
}
dominant_structure_for_cluster <- function(df, rename){
  my_cluster <- function(x){
    result_of_cluster <- tsclust(
      x,
      k = 6,                 # запрашиваемое число кластеров
      type = "hierarchical", # тип кластеризации
      distance = "dtw",      # мера расстояния - динамическая трансформации временной шкалы (Dynamic Time Warping)
      seed = 42,
      control = 
        hierarchical_control(method = "ward.D2"), # метод агломерации
      args = 
        tsclust_args(dist = list(window.size = 7)) # размер окна Сакэ-Чиба
    )
    return(result_of_cluster)
  }#запустить кластеризацию
  plot_of_cluster <- function(x, names){
    return(plot(x, xlab = "", sub = "", main = names))
  }
  df <- df%>%
    mutate_at(vars(c(2:ncol(df))), funs(.*100/df[,ncol(df)]))%>%
    mutate_if(is.numeric, round, digits=1)%>%
    select(-Community)%>%
    rename_with(~paste0(rename, .x), .cols=-1)
  df <- as.list(df[2:ncol(df)])
  cluster_of_number <- my_cluster(df) 
  dendrogram_of_numbers <- plot_of_cluster(cluster_of_number, "Кластеризация по структуре доминирования")
  return(dendrogram_of_numbers)
  
} #Кластеризация по  структуре доминирования
stationar <- function(x){
  afd <- function(x){
    a <- adf.test(x)
    if (a$p.value < 0.05){
      return(a)
    }
  }
  return(apply(x[2:ncol(x)], 2, afd))
} #Проверка ряда на стационарность тестом Дики-Фулера

################################################################################
############# Объединенные ряды XX и XXI века   ################################

wavelet_cluster(LB, "LB_")

?analyze.wavelet





################################################################################
##########               XX век            #####################################

# Результаты вейвлет-преобразования и кластеризации
result_LB_XX <- wavelet_cluster(LB_XX, "LB_")
result_RB_XX <- wavelet_cluster(RB_XX, "RB_")
result_both_bank_XX <- wavelet_cluster(XX_both_bank, "All_") #суммарная численность


#автокорреляция
test_autocorrelation(LB_XX[2:ncol(LB_XX)])
test_autocorrelation(RB_XX[2:ncol(RB_XX)])
test_autocorrelation(XX_both_bank[2:ncol(XX_both_bank)])


################################################################################
##########               XXI век            ####################################


#автокорреляция - строгих циклов в 21 веке не существует
test_autocorrelation(LB_XXI[2:ncol(LB_XXI)])
test_autocorrelation(RB_XXI[2:ncol(RB_XXI)])
test_autocorrelation(XXI_both_bank[2:ncol(XXI_both_bank)])



# проверка на стационарность
# стационарный временной ряд не имеет тренда
# нужно учитывать, что если ряд не стационарный, то для анализа временных рядов
# его нужно сначала лишить тренда

stationar(LB_XXI) #стационарные ряды для левого берега: community, sic_bet,
# Cle_rfc, Tal_alt, Sor_mss
stationar(RB_XXI)#стационарные ряды для правого берега: community, sic_bet,
# cle_rfc, Sor_iso, sor_ara
stationar(XXI_both_bank)#стационарные ряды для правого берега: community, sic_bet,
# cle_rut, Sor_iso

# Но у многих видов начинает появляться что-то похожее на цикл! - см вейвлет
# Но что это - покажут дальнейшие учеты





# Результаты вейвлет-преобразования и кластеризации
result_LB_XXI <- wavelet_cluster(LB_XXI, "LB_")
result_RB_XXI <- wavelet_cluster(RB_XXI, "RB_")
result_both_bank_XXI <- wavelet_cluster(XXI_both_bank, "All_") #суммарная численность


#расчет доли вида в сообществе - индекс доминирования
XXI_both_bank_dominant <- dominant_structure_for_graph(XXI_both_bank)
XXI_LB_dominant <- dominant_structure_for_graph(LB_XXI)
XXI_RB_dominant <- dominant_structure_for_graph(RB_XXI)


#Попробуем как-то классифицировать виды по динамике структуры доминирования
dominant_structure_for_cluster(LB_XXI, 'LB_')
dominant_structure_for_cluster(RB_XXI, 'RB_')
dominant_structure_for_cluster(XXI_both_bank, 'ALL_')

#Графики численности и структуры доминирования
my_graph <- function(df){
  graph <- function(df, col){
    return(ggplot(df)+
             geom_line(aes(date, col))+
             geom_point(aes(date, col))+
             labs(x = "Год" , y = "Численность (экземляров на 100 ц-с", title = col)
    )
  }
  return(map(df[2:ncol(df)], graph, df = df))
}



my_graph(XX_both_bank)

XX_community <- XX_both_bank%>%
  filter(date >=1982 & date<=1990)%>%
  mutate(date = as.numeric(date))
ggplot(XX_community, aes(date, Community))+
  geom_line()+
  geom_point()+
  labs(x = "Год" , y = "Численность (экземляров на 100 ц-с")+
  theme_classic()

ggplot(XXI_both_bank, aes(date, Community))+
  geom_line()+
  geom_point()+
  labs(x = "Год" , y = "Численность (экземляров на 100 ц-с")+
  theme_classic()


#на левом берегу динамика этих видов
# практически точно следует за численностью сообщества





# Настроить шкалы графика
# считать индексы доминирования
# индексы разнообразия 
# делать временные ряды и их тоже анализировать
# Как связана доля вида в сообществе и численность сообщества? (в зависимости от вида)
# Каким-то образом классифицировать виды (по численности, по )


# Что дальше
# Читать про интерпретацию вейвлет преобразования - что если период дробный,
# а разрешение метода - строго целое число? (как в нашем случае, разрешение
# - год, а период бывает промежуточным)
# Тонко настроить график, выбрать функции и тд

# Кластеризация
# Число кластеров
# Достоверность подразделения на клады - такой не бывает, но может быть матрица расстояний?
# Объединить в общий список суммарную численность и отдельно по берегам
# (добавить индексы к названиям столбцов)
# в итоге получим группы видов со сходной периодичностью
#отдельно сделать кластеризацию исходных рядов численности 
# (можно логарифмировать для снижения дисперсии)
#в итоге получим похожесть рядов друг на друга
#при этом в мере расстояния нужно будет указать евклидово расстояние
# тогда по сути будет оценена синхронность, ведь даже при небольшом смещении врем рядов
# Друг относительно друга будет другой кластер - Это только предположение

#визуализация 
# общая картинка со спектрами для всех видов (делим по берегам)
#общая картинка по ходу численности для всех видов по берегам


?tsclust
#Попробовать поискать точки излома https://ranalytics.github.io/tsa-with-r/ch-structural-changes.html


