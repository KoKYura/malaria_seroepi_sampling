library(dplyr)
library(tidyverse)
library(shiny)
library(shinydashboard)
library(DT)
library(JuliaCall)

showtext::showtext_auto(enable = TRUE)
options(dplyr.summarise.inform=F)
julia_setup()
julia_command(paste0("include(\"", getwd(), "/functions.jl\")"))

# Get possible age allocations for each scenario
df_age_comb_s1 <- read.csv("./data/df_age_comb_s1.csv")
df_age_comb_s2 <- read.csv("./data/df_age_comb_s2.csv")

# Functions in R ----

col_q1 <- c("white",'#ffeda0','#feb24c','#fd8d3c','#fc4e2a','#bd0026')
col_q5 <- c("white",'#deebf7','#9ecae1','#4292c6','#08519c','#08306b')

get_ci <- function(df1, df2, target_est_, SCR_, title_){
  df_summary_q1 <- df1 %>% filter(score == "Q1") %>%  dplyr::select(i,score)
  df_summary_q5 <- df1 %>% filter(score == "Q5") %>%  dplyr::select(i,score)
  
  df_est_q1 <- left_join(df2, df_summary_q1, by = "i")  %>%  filter(!is.na(score))
  df_est_q5 <- left_join(df2, df_summary_q5, by = "i")  %>%  filter(!is.na(score))
  
  df_ci <- rbind(df_est_q1,df_est_q5) %>% 
    group_by(score) %>% 
    rename(target_est = target_est_) %>% 
    mutate(est_l_med = median(target_est),
           est_l_upper = quantile(target_est, 0.975),
           est_l_lower = quantile(target_est, 0.025)) %>% 
    distinct(score, .keep_all = T) 
  
  g_ci <- ggplot(df_ci, aes(x = as.factor(score), y = est_l_med, color = score)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(ymin = est_l_lower, ymax = est_l_upper), 
                  position = position_dodge(width = 0.5), 
                  size = 1,
                  width = 0.4)+
    theme(panel.grid.major.y = element_blank(),
          axis.text = element_text(size = 10))+
    labs(y = "Estimated SCR", x = "") +
    theme_bw()+
    ggtitle(title_)+
    theme(plot.title = element_text(size = 14, face = "bold"),
          panel.grid.major.y = element_blank(),
          axis.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12, face = "bold"),
          strip.background = element_rect(fill = "white", colour = "black")
    )+
    geom_hline(yintercept = SCR_, linetype = "dashed")
  g_ci
}
get_heat <- function(df, df_age_comb_, score_, tau_, col_, title_){
  
  
  age1_c <- paste0("[1,",tau_,"]")
  age2_c <- paste0("[",1+tau_,",",10+tau_,"]")
  age3_c <- paste0("[",11+tau_,",",20+tau_,"]")
  age4_c <- paste0("[",21+tau_,",",30+tau_,"]")
  age5_c <- paste0("[",31+tau_,",",40+tau_,"]")
  
  df_score <- df %>% 
    filter(score == score_) %>% 
    left_join(df_age_comb_, by = "i") %>% 
    arrange(age1,age2,age3,age4) %>% 
    mutate(order = factor(row_number(), levels = row_number())) %>% 
    dplyr::select(order, age1, age2, age3, age4, age5) %>% 
    pivot_longer(cols = starts_with("age"),
                 names_to  = "age_group",
                 values_to = "prop") %>% 
    mutate(age_group = factor(age_group, levels = paste0("age", 1:5))) %>% 
    mutate(prop2 = case_when(prop < 0.1 ~ "0%",
                             prop < 0.3 ~ "20%",
                             prop < 0.5 ~ "40%",
                             prop < 0.7 ~ "60%",
                             prop < 0.9 ~ "80%",
                             prop >= 0.9 ~ "100%"),
           age_group2 = case_when(age_group == "age1" ~ age1_c,
                                  age_group == "age2" ~ age2_c,
                                  age_group == "age3" ~ age3_c,
                                  age_group == "age4" ~ age4_c,
                                  age_group == "age5" ~ age5_c))
  
  df_score$prop2 <- fct_relevel(df_score$prop2, c("0%", "20%", "40%", "60%", "80%", "100%"))
  df_score$age_group2 <- fct_relevel(df_score$age_group2, c(age1_c, age2_c, age3_c, age4_c, age5_c))
  
  ggplot(df_score, aes(x = age_group2, y = order, fill = prop2)) +
    scale_fill_manual(values = col_)+
    geom_tile(colour = "grey60") +
    labs(x = "Age group", y = "") +
    theme_bw()+
    ggtitle(title_)+
    theme(plot.title = element_text(size = 14, face = "bold"),
          panel.grid.major.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          strip.text = element_text(size = 12, face = "bold"),
          strip.background = element_rect(fill = "white", colour = "black")
    )
  
}

# UI ----
ui <- dashboardPage(
  dashboardHeader(title = "Optimal Age-Sampling Structures for Malaria Seroprevalence",
                  titleWidth = 500),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("Stable SCR", tabName = "page1"),
                menuItem("Change in SCR", tabName = "page2")
    )
  ),
  dashboardBody(
    fluidRow(
        tags$style(HTML("
        .main-header .logo {
        font-weight: bold !important;
        background-color: #252525 !important;
        }
        .main-header .navbar {
        background-color: #4d4d4d !important;
        }
        .main-sidebar {
        background-color: #222 !important; 
        }
        .content-wrapper {
        padding-left: 30px !important;
        padding-right: 30px !important;
        background-color: #878787 !important; 
      　}
        .box-title {
        font-weight: bold; 
        }
        input[type='number']::-webkit-outer-spin-button,
        input[type='number']::-webkit-inner-spin-button {
          -webkit-appearance: none;
          margin: 0;
        }
      ")),
    tabItems(
      tabItem(tabName = "page1",
              fluidRow(
                column(6,
                       wellPanel(
                         numericInput("n_rec", h4("Number of samples"), value = NA, min = 250, max = 5000),
                         numericInput("SCR", h4("SCR"), value = NA, min = 0, max = 1),
                         numericInput("SRR", h4("SRR"), value = NA, min = 0, max = 1),
                         selectInput("param", h4("Estimating Parameter"), choices = c("SCR only", "SCR and SRR"), selected = NULL)
                       )
                     )
                   ),
          fluidRow(
            column(12, div(
              style = "margin-bottom: 15px;",
              actionButton("run1", "Run", 
                           style = "background-color: #252525; color: white;")
            ))
          ),
          fluidRow(
            box(
              title = "",
              width = 12,
              plotOutput("plot1")
            )
          )
        ),
      tabItem(tabName = "page2",
              fluidRow(
                column(6,
                       wellPanel(
                         numericInput("n_rec_s2", h4("Number of samples"), value = NA, min = 250, max = 5000),
                         numericInput("SCR1_s2", h4("SCR before the change point"), value = NA, min = 0, max = 1),
                         numericInput("SCR2_s2", h4("SCR after the change point"), value = NA, min = 0, max = 1),
                         numericInput("SRR_s2", h4("SRR"), value = NA, min = 0, max = 1),
                         sliderInput("tau", h4("Change point (Years)"), min = 1, max = 20, value = 5, step = 1),
                         selectInput("param_s2", h4("Estimating Parameters"), choices = c("SCRs only", "SCRs and SRR"), selected = NULL)
                       )
                )
              ),
              fluidRow(
                column(12, div(
                  style = "margin-bottom: 15px;",
                  actionButton("run2", "Run", 
                               style = "background-color: #252525; color: white;")
                ))
              ),
              fluidRow(
                box(
                  title = "",
                  width = 12,
                  plotOutput("plot2")
                )
              )
        )
      )
    )
  )
)

# Server ----
server <- function(input, output, session) {
  
  # Scenario 1
  
  input_ready_s1 <- reactiveVal(FALSE)
  observe({
    if (!is.null(input$n_rec) && !is.null(input$SCR) && !is.null(input$SRR)) {
      input_ready_s1(TRUE)
    }
  })
  
  observeEvent(input$run1, {
    req(input_ready_s1())
    output$plot1 <- renderPlot({
      
      n_rec <- input$n_rec
      SCR <- input$SCR
      SRR <- input$SRR
      param <- input$param
        
      withProgress(message = "Calculating", value = 0, {
        
      df_est <- tibble()
      for (i in 1:nrow(df_age_comb_s1)){
        Sys.sleep(0.1)
        setProgress(value = i/nrow(df_age_comb_s1),
                    detail = sprintf("%d/%d", i, nrow(df_age_comb_s1)))
        
        p_age1 = df_age_comb_s1[i,][[2]]; p_age2 = df_age_comb_s1[i,][[3]]; p_age3 = df_age_comb_s1[i,][[4]]; p_age4 = df_age_comb_s1[i,][[5]]
        
        if (param == "SCR only"){
          res <- lapply(1:1000, function(x) {
            julia_call("get_res_s1_par1", as.integer(n_rec), p_age1, p_age2, p_age3, p_age4, SCR, SRR)
          })
           
        }else{
          res <- lapply(1:1000, function(x) {
            julia_call("get_res_s1_par2", as.integer(n_rec), p_age1, p_age2, p_age3, p_age4, SCR, SRR)
          }) 
        }
        
        df_est_ <- do.call(rbind, res) %>% 
          mutate(i = i)
        
        df_est <- rbind(df_est, df_est_)
       }
      })
      
      df_summary <- df_est %>% 
        group_by(i) %>% 
        mutate(prec_lambda = (quantile(est_l,0.975) - quantile(est_l,0.025))/2) %>% 
        mutate(rel_prec = prec_lambda/SCR) %>% 
        ungroup() %>% 
        distinct(i, .keep_all = T) %>% 
        mutate(score = case_when(rel_prec < quantile(rel_prec, 0.20) ~ "Q1",
                                 rel_prec < quantile(rel_prec, 0.40) ~ "Q2",
                                 rel_prec < quantile(rel_prec, 0.60) ~ "Q3",
                                 rel_prec < quantile(rel_prec, 0.80) ~ "Q4",
                                 rel_prec >= quantile(rel_prec, 0.80) ~ "Q5"
        ))
      
      g_ci_s1 <- get_ci(df_summary, df_est, "est_l", SCR, "95% CI for the SCR")
      g_heat_Q1_s1 <- get_heat(df_summary, df_age_comb_s1, "Q1", 10, col_q1, "Optimal age structures")
      g_heat_Q5_s1 <- get_heat(df_summary, df_age_comb_s1, "Q5", 10, col_q5, "Least precise age structures")
  
      cowplot::plot_grid(g_ci_s1, NA, g_heat_Q1_s1, g_heat_Q5_s1, nrow = 2)
    })
  })
      
  # Scenario 2
  
  input_ready_s2 <- reactiveVal(FALSE)
  
  observe({
    if (!is.null(input$n_rec_s2) && !is.null(input$SCR1_s2) && !is.null(input$SCR2_s2) && 
        !is.null(input$SRR_s2) && !is.null(input$tau) && !is.null(input$param_s2)) {
      input_ready_s2(TRUE)
    }
  })
  
  observeEvent(input$run2, {
    req(input_ready_s2())
    output$plot2 <- renderPlot({
      n_rec_s2 <- input$n_rec_s2
      SCR1_s2 <- input$SCR1_s2
      SCR2_s2 <- input$SCR2_s2
      SRR_s2 <- input$SRR_s2
      tau <- input$tau
      param_s2 <- input$param_s2
      
      withProgress(message = "Calculating", value = 0, {
        
        df_est_s2 <- tibble()
        for (i in 1:nrow(df_age_comb_s2)){
          Sys.sleep(0.1)
          setProgress(value = i/nrow(df_age_comb_s2),
                      detail = sprintf("%d/%d", i, nrow(df_age_comb_s2)))
          
          p_age1 = df_age_comb_s2[i,][[2]]; p_age2 = df_age_comb_s2[i,][[3]]; p_age3 = df_age_comb_s2[i,][[4]]; p_age4 = df_age_comb_s2[i,][[5]]
          
          if (param_s2 == "SCRs only"){
            res <- lapply(1:1000, function(x) {
              julia_call("get_res_s2_par2", as.integer(n_rec_s2), p_age1, p_age2, p_age3, p_age4, SCR1_s2, SCR2_s2, SRR_s2, tau)
            })
            
          }else{
            res <- lapply(1:1000, function(x) {
              julia_call("get_res_s2_par3", as.integer(n_rec_s2), p_age1, p_age2, p_age3, p_age4, SCR1_s2, SCR2_s2, SRR_s2, tau)
            }) 
          }
          
          df_est_ <- do.call(rbind, res) %>% 
            mutate(i = i)
          
          df_est_s2 <- rbind(df_est_s2, df_est_)
        }
      })
      
      df_summary_s2 <- df_est_s2 %>% 
        group_by(i) %>% 
        mutate(prec_l1 = (quantile(est_l1,0.975) - quantile(est_l1,0.025))/2,
               prec_l2 = (quantile(est_l2,0.975) - quantile(est_l2,0.025))/2) %>% 
        mutate(joint_prec = prec_l1/SCR1_s2 + prec_l2/SCR2_s2) %>% 
        ungroup() %>% 
        distinct(i, .keep_all = T) %>% 
        mutate(score = case_when(joint_prec < quantile(joint_prec, 0.20) ~ "Q1",
                                 joint_prec < quantile(joint_prec, 0.40) ~ "Q2",
                                 joint_prec < quantile(joint_prec, 0.60) ~ "Q3",
                                 joint_prec < quantile(joint_prec, 0.80) ~ "Q4",
                                 joint_prec >= quantile(joint_prec, 0.80) ~ "Q5"))
      
      g_ci_l1_s2 <- get_ci(df_summary_s2, df_est_s2, "est_l1", SCR1_s2, "95% CI for the SCR before the change")
      g_ci_l2_s2 <- get_ci(df_summary_s2, df_est_s2, "est_l2", SCR2_s2, "95% CI for the SCR after the change")
      g_heat_Q1_s2 <- get_heat(df_summary_s2, df_age_comb_s2, "Q1", tau, col_q1, "Optimal age structures")
      g_heat_Q5_s2 <- get_heat(df_summary_s2, df_age_comb_s2, "Q5", tau, col_q5, "Least precise age structures")
      
      cowplot::plot_grid(g_ci_l1_s2, g_ci_l2_s2, g_heat_Q1_s2, g_heat_Q5_s2, nrow = 2)
    })
  })
}

# Run app ----
shinyApp(ui,server)





