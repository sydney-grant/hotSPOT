library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(dashboardthemes)
library(DT)
library(plotly)


customTheme <- shinyDashboardThemeDIY(
    
    ### general
    appFontFamily = "Arial"
    ,appFontColor = "rgb(50,50,50)"
    ,primaryFontColor = "rgb(50,50,50)"
    ,infoFontColor = "rgb(0,0,0)"
    ,successFontColor = "rgb(0,0,0)"
    ,warningFontColor = "rgb(0,0,0)"
    ,dangerFontColor = "rgb(0,0,0)"
    ,bodyBackColor = "rgb(250,250,250)"
    
    ### header
    ,logoBackColor = "rgb(255,191,122)"
    ,headerButtonBackColor = "rgb(255,191,122)"
    ,headerButtonIconColor = "rgb(50,50,50)"
    ,headerButtonBackColorHover = "rgb(255,191,122)"
    ,headerButtonIconColorHover = "rgb(50,50,50)"
    
    ,headerBackColor = cssGradientThreeColors(
        direction = "right"
        ,colorStart = "rgba(255,191,122,1)"
        ,colorMiddle = "rgba(255,242,122,1)"
        ,colorEnd = "rgba(252, 242, 141,1)"
        ,colorStartPos = 0
        ,colorMiddlePos = 90
        ,colorEndPos = 100
    )
    ,headerBoxShadowColor = "#aaaaaa"
    ,headerBoxShadowSize = "2px 2px 2px"
    
    ### sidebar
    ,sidebarBackColor = cssGradientThreeColors(
        direction = "down"
        ,colorStart = "rgb(255,191,122)"
        ,colorMiddle = "rgb(255,242,122)"
        ,colorEnd = "rgb(255,138,122)"
        ,colorStartPos = 0
        ,colorMiddlePos = 50
        ,colorEndPos = 100
    )
    ,sidebarPadding = 0
    
    ,sidebarMenuBackColor = "transparent"
    ,sidebarMenuPadding = 0
    ,sidebarMenuBorderRadius = 0
    
    ,sidebarShadowRadius = "3px 5px 5px"
    ,sidebarShadowColor = "#aaaaaa"
    
    ,sidebarUserTextColor = "rgb(255,255,255)"
    
    ,sidebarSearchBackColor = "rgb(55,72,80)"
    ,sidebarSearchIconColor = "rgb(153,153,153)"
    ,sidebarSearchBorderColor = "rgb(55,72,80)"
    
    ,sidebarTabTextColor = "rgb(50,50,50)"
    ,sidebarTabTextSize = 13
    ,sidebarTabBorderStyle = "none none solid none"
    ,sidebarTabBorderColor = "rgb(35,106,135)"
    ,sidebarTabBorderWidth = 1
    
    ,sidebarTabBackColorSelected = cssGradientThreeColors(
        direction = "right"
        ,colorStart = "rgba(255,138,122,1)"
        ,colorMiddle = "rgba(252,110,90,1)"
        ,colorEnd = "rgba(212,90,70,1)"
        ,colorStartPos = 0
        ,colorMiddlePos = 30
        ,colorEndPos = 100
    )
    ,sidebarTabTextColorSelected = "rgb(250,250,250)"
    ,sidebarTabRadiusSelected = "0px 20px 20px 0px"
    
    ,sidebarTabBackColorHover = cssGradientThreeColors(
        direction = "right"
        ,colorStart = "rgba(255,242,122,1)"
        ,colorMiddle = "rgba(252, 242, 141,1)"
        ,colorEnd = "rgba(255, 248, 179,1)"
        ,colorStartPos = 0
        ,colorMiddlePos = 30
        ,colorEndPos = 100
    )
    ,sidebarTabTextColorHover = "rgb(50,50,50)"
    ,sidebarTabBorderStyleHover = "none none solid none"
    ,sidebarTabBorderColorHover = "rgb(75,126,151)"
    ,sidebarTabBorderWidthHover = 1
    ,sidebarTabRadiusHover = "0px 20px 20px 0px"
    
    ### boxes
    ,boxBackColor = "rgb(255,255,255)"
    ,boxBorderRadius = 5
    ,boxShadowSize = "0px 1px 1px"
    ,boxShadowColor = "rgba(0,0,0,.1)"
    ,boxTitleSize = 16
    ,boxDefaultColor = "rgb(212,90,70)"
    ,boxPrimaryColor = "rgba(212,90,70,1)"
    ,boxInfoColor = "rgb(255, 252, 166)"
    ,boxSuccessColor = cssGradientThreeColors(
        direction = "right"
        ,colorStart = "rgba(255,138,122,1)"
        ,colorMiddle = "rgba(252,110,90,1)"
        ,colorEnd = "rgba(212,90,70,1)"
        ,colorStartPos = 0
        ,colorMiddlePos = 30
        ,colorEndPos = 100)
    ,boxWarningColor = "rgb(244,156,104)"
    ,boxDangerColor = "rgb(255,88,55)"
    
    ,tabBoxTabColor = "rgb(212,90,70,1)"
    ,tabBoxTabTextSize = 14
    ,tabBoxTabTextColor = "rgb(0,0,0)"
    ,tabBoxTabTextColorSelected = "rgb(0,0,0)"
    ,tabBoxBackColor = "rgb(212,90,70,1)"
    ,tabBoxHighlightColor = "rgba(212,90,70,1,1)"
    ,tabBoxBorderRadius = 5
    
    ### inputs
    ,buttonBackColor = cssGradientThreeColors(
        direction = "right"
        ,colorStart = "rgba(255,138,122,1)"
        ,colorMiddle = "rgba(252,110,90,1)"
        ,colorEnd = "rgba(212,90,70,1)"
        ,colorStartPos = 0
        ,colorMiddlePos = 30
        ,colorEndPos = 100)
    ,buttonTextColor = "rgb(250,250,250)"
    ,buttonBorderColor = "rgb(200,200,200)"
    ,buttonBorderRadius = 5
    
    ,buttonBackColorHover = cssGradientThreeColors(
        direction = "right"
        ,colorStart = "rgba(255,242,122,1)"
        ,colorMiddle = "rgba(252, 242, 141,1)"
        ,colorEnd = "rgba(255, 248, 179,1)"
        ,colorStartPos = 0
        ,colorMiddlePos = 30
        ,colorEndPos = 100
    )
    ,buttonTextColorHover = "rgb(50,50,50)"
    ,buttonBorderColorHover = "rgb(200,200,200)"
    
    ,textboxBackColor = "rgb(255,255,255)"
    ,textboxBorderColor = "rgb(200,200,200)"
    ,textboxBorderRadius = 5
    ,textboxBackColorSelect = "rgb(245,245,245)"
    ,textboxBorderColorSelect = "rgb(200,200,200)"
    
    ### tables
    ,tableBackColor = "rgb(255,255,255)"
    ,tableBorderColor = "rgb(240,240,240)"
    ,tableBorderTopSize = 1
    ,tableBorderRowSize = 1
    
)

customLogo <- shinyDashboardLogoDIY(
    boldText = "hotSPOT"
    ,mainText = "(Sequencing. Panel. Optimization. Tool.)"
    ,textSize = 20
    ,badgeText = "v1.1"
    ,badgeTextColor = "black"
    ,badgeTextSize = 1
    ,badgeBackColor = "#ffbf7a"
    ,badgeBorderRadius = 3
    
)

ui <- dashboardPage(
    
    dashboardHeader(title = customLogo, 
                    titleWidth = 500),
    dashboardSidebar(sidebarMenu(
        menuItem("Forward Selection Method", tabName = "forward"),
        menuItem("Comprehensive Selection Method", tabName = "exhaustive")
    )),
    dashboardBody(customTheme,
                  tabItems(
        tabItem(tabName = "forward", h4("Targeted Sequencing Panel Design: Forward Selection Method"),
                fluidRow(
                    box(title = "Input Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                fileInput(inputId = "data1",
                          label = "Data Set:",
                          accept = ".csv",
                          width = NULL),
                numericInput(inputId = "panel_length",
                             label = "Length of Sequencing Panel (#bp):",
                             10000),
                numericInput(inputId = "amp_length",
                             label = "Length of Amplicons (#bp):",
                             125),
                radioButtons("radio", label = "Include Genes:",
                             choices = list("Yes" = 1, "No" = 2), 
                             selected = 1),
                actionButton("goButton", "Submit"), 
                downloadButton("myFile", label = "Download Panel"), height = 425),
                box(title = "Panel Capture Efficacy", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    plotlyOutput("myPlot", height = 363) %>% withSpinner(color = "black"))),
                fluidRow(box(title = "Output Summary:", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                             textOutput("text"))),
                dataTableOutput("myData") %>% withSpinner(color = "black")),
        tabItem(tabName = "exhaustive", h4("Targeted Sequencing Panel Design: Comprehensive Selection Method"),
                fluidRow(
                    box(title = "Input Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                fileInput(inputId = "data2",
                          label = "Data Set:",
                          accept = ".csv",
                          width = NULL),
                numericInput(inputId = "panel_length2",
                             label = "Length of Sequencing Panel (#bp):",
                             10000),
                numericInput(inputId = "amp_length2",
                             label = "Length of Amplicons (#bp):",
                             125),
                sliderInput(inputId = "num_amps",
                             label = "Size of Hotspots(#amplicons):",
                             min = 1, max = 10, value = 3),
                radioButtons("radio", label = "Include Genes:",
                             choices = list("Yes" = 1, "No" = 2), 
                             selected = 1),
                actionButton("goButton2", "Submit"), 
                downloadButton("myFile2", label = "Download Panel"), height = 525),
                box(title = "Panel Capture Efficacy", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    plotlyOutput("myPlot2", height = 464) %>% withSpinner(color = "black"))),
                fluidRow(box(title = "Output Summary:", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                             textOutput("text2"))),
                dataTableOutput("myData2") %>% withSpinner(color = "black"))
    )
)
)


server <- function(input, output) {

        myPlot = reactiveVal()
        myData = reactiveVal()
        myFile = reactiveVal()
        text = reactiveVal()
        
        observeEvent(input$goButton, {
        
        data1 <- req(input$data1)
        data <- read.csv(data1$datapath)
        
        panel_length = req(input$panel_length) 
        amp_length = req(input$amp_length)
        
        ######### amplicon based binning algorithm
        
        library(dplyr)
        library(plyr)
        library(hash)
        library(rlist)
        library(R.utils)
        library(ggplot2)
        library(plotly)
        
        pos <- data$pos
        pos_freq <- data.frame(ftable(pos)) #make frequency table for each position
        
        data <- merge(data, pos_freq, by = "pos") #merge frequency table with original df
        data <- unique(data) #keep only unique positions
        
        
        amplicon_finder <- function(data, chr){
          pos <- data$pos
          
          keys <- make.keys(data$pos)
          dict <- hash(keys = keys, values = data$Freq) #create dictionary for mutation locations and mutation frequency
          
          make_bins <- sort(unique(pos))
          bins <- data.frame()
          while (length(make_bins) > 1){
            temp_bins <- make_bins
            i <- 1
            start <- make_bins[1]
            p_dif <- 1
            go = "yes"
            while (go == "yes"){
              p1 <- temp_bins[i]
              p_next <- temp_bins[i + 1]
              i <- i + 1
              if (i <= length(make_bins)){
                p_dif <- p_next - p1
              }
              if (i > length(make_bins)){
                p_dif = amp_length + 100
                go = "no"
              }
              if (p_dif > amp_length){
                p_next <- p1
                go = "no"
              }
            }
            end <- p_next
            vec <- as.list(start:end)
            keep <- setdiff(make_bins, vec)
            make_bins <- keep
            chr = chr
            row <- data.frame("lowerbound" = start, "upperbound" = end, "chromosome" = chr)
            bins <- rbind(bins, row)
          }
          if (length(make_bins) == 1){
            row <- data.frame("lowerbound" = make_bins, "upperbound" = make_bins, "chromosome" = chr)
            bins <- rbind(bins, row)
          }
          
          mut_locations <- data.frame()
          possible_bins <- data.frame()
          
          for (i in 1:nrow(bins)){
            up_mut <- bins$upperbound[[i]]
            low_mut <- bins$lowerbound[[i]]
            vec <- up_mut:low_mut
            dif <- length(vec) - amp_length
            if (dif <= 0){ #############FIX HERE#######################
              up <- bins$upperbound[[i]]
              low <- bins$lowerbound[[i]]
              vec <- up:low
              mutations <- intersect(vec, pos) #find all mutations in total bin region
              mutations_keys <- make.keys(mutations)
              count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
              mut_all <- c()
              for (m in 1:length(mutations)){
                mut_all <- c(mut_all, rep(mutations[[m]], as.vector(values(dict, keys = mutations_keys))[[m]]))
              }
              weighted_mid <- round(mean(mut_all))
              weighted_mid_max <- (ceiling(mean(vec)) + (amp_length/2))
              weighted_mid_min <- (floor(mean(vec)) - (amp_length/2))
              if (weighted_mid < weighted_mid_min){final_mid <- weighted_mid_min}
              if (weighted_mid > weighted_mid_max){final_mid <- weighted_mid_max}
              if (weighted_mid > weighted_mid_min & weighted_mid < weighted_mid_max){final_mid <- weighted_mid}
              up <- final_mid + ceiling((amp_length - 1) / 2)
              low <- final_mid - floor((amp_length - 1) / 2)
              vec <- up:low
              mutations <- intersect(vec, pos) #find all mutations in total bin region
              mutations_keys <- make.keys(mutations)
              count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
              row <- data.frame("lowerbound" = low, "upperbound" = up, 
                                "chromosome" = bins$chromosome[[i]], "count" = count,"id" = "x")
              possible_bins <- rbind(possible_bins, row)
            }
            else{
              up_mut <- bins$upperbound[[i]]
              low_mut <- bins$lowerbound[[i]]
              vec <- up_mut:low_mut
              mut_list <- intersect(vec, pos)
              id = paste(as.character(bins$chromosome[[i]]), as.character(i), sep = "-")
              b_list <- list()
              for (mut in mut_list){
                bin1 <- (mut - amp_length):(mut - 1)
                b_list[[1]] <- bin1
                bin2 <- (mut - (amp_length - 1)):(mut)
                b_list[[2]] <- bin2
                bin3 <- (mut + amp_length):(mut + 1)
                b_list[[3]] <- bin3
                bin4 <- (mut + (amp_length - 1)):(mut)
                b_list[[4]] <- bin4
                bin5 <- (mut - ceiling(amp_length / 2)):(mut - floor(amp_length / 2))
                b_list[[5]] <- bin5
                
                for (b in b_list){
                  low <- min(b)
                  up <- max(b)
                  mutations <- intersect(b, pos)
                  if (length(mutations > 0)){
                    mutations_keys <- make.keys(mutations)
                    count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                    row <- data.frame("lowerbound" = low, "upperbound" = up, 
                                      "chromosome" = bins$chromosome[[i]], "count" = count,"id" = id)
                    possible_bins <- rbind(possible_bins, row)
                  }
                }
              }
              
            }
          }
          return(possible_bins)
        }
        
        
        
        
        ##################################################################################
        
        
        chrom_list <- unique(as.list(data$chr)) # list of all chromosomes found in this dataset
        all_pos_bins <- data.frame()
        
        for (chrom in chrom_list) {
          chromosome <- subset(data, chr == chrom)
          chromsort <- chromosome[order(chromosome$pos), ] #order df from lowest to highest position
          chrom_bins <- amplicon_finder(chromsort, chrom) #run through above function
          all_pos_bins <- rbind(all_pos_bins, chrom_bins) #add all bins together from all chromosomes
        }
        
        ###################################################################################
        ordered_bin <- all_pos_bins[order(-all_pos_bins$count), ]
        #ordered_bin <- ordered_bin[order(nrow(ordered_bin):1),]
        
        big_spots <- ordered_bin[!ordered_bin$id == "x", ]
        unique_spots <- unique(as.list(big_spots$id))
        
        keys <- make.keys(data$pos)
        dict <- hash(keys = keys, values = data$Freq)
        possible_bins <- data.frame()
        for (spot in unique_spots){
          spot_sub <- subset(big_spots, id == spot)
          big_bin_all <- list()
          while (nrow(spot_sub) > 0){
            row <- data.frame()
            new_sub <- data.frame()
            big_spot <- spot_sub[1,]
            big_bin <- big_spot$upperbound:big_spot$lowerbound
            big_bin_all <- c(big_bin_all, big_bin)
            possible_bins <- rbind(possible_bins, big_spot)
            if (nrow(spot_sub) > 1){
              for (k in 2:nrow(spot_sub)){
                row <- data.frame()
                test_bin <- spot_sub$upperbound[k]:spot_sub$lowerbound[k]
                overlap <- intersect(big_bin_all, test_bin)
                if (length(overlap > 0)){
                  unique_bin <- setdiff(test_bin, overlap)
                  mutations <- intersect(unique_bin, pos) #find all mutations in total bin region
                  if (length(mutations) > 0){
                    mutations_keys <- make.keys(mutations)
                    new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                    row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin), 
                                      "chromosome" = spot_sub$chromosome[1], "count" = new_count,"id" = spot_sub$id[1])
                  }
                }
                if (length(overlap) <= 0) {row <- spot_sub[k,]}
                new_sub <- rbind(new_sub, row)
              }
            }
            spot_sub <- new_sub
            if (nrow(spot_sub > 0)){spot_sub <- spot_sub[order(-spot_sub$count), ]}
          }
        }
        
        
        final_bin <- rbind(subset(ordered_bin, id == "x"), possible_bins)
        final_bin <- final_bin[order(-final_bin$count), ]
        
        # reorder to favor amplicons with mutations in center
        
        cutoff_point <- final_bin$count[[round(panel_length / amp_length)]]
        
        subset_final_bin <- subset(final_bin, count == cutoff_point)
        dif_list <- list()
        for (i in nrow(subset_final_bin)){
          whole_region <- subset_final_bin$lowerbound[[i]]:subset_final_bin$upperbound[[i]]
          bin_mutations <- intersect(pos, whole_region)
          mid_diff <- abs((max(whole_region) - max(bin_mutations)) - (min(bin_mutations) - min(whole_region)))
          dif_list <- c(dif_list, mid_diff)
        }
        subset_final_bin$dif <- unlist(dif_list)
        subset_final_bin <- subset_final_bin[order(-subset_final_bin$dif), ]
        subset_final_bin <- subset_final_bin[,-ncol(subset_final_bin)]
        
        df_range_min <- min(which(final_bin$count == cutoff_point))
        df_range_max <- max(which(final_bin$count == cutoff_point))
        
        final_bin[df_range_min:df_range_max,] <- subset_final_bin[1:nrow(subset_final_bin),]
        
        ##############################################################
        
        
        bin_len_list <- list()
        for (y in 1:nrow(final_bin)) {      #calculating bin length using upper and lower bound positions
          upper = final_bin$upperbound[[y]]
          lower = final_bin$lowerbound[[y]]
          bin_len <- upper - lower + 1# add 1 to include both upper and lowerbound plus all bp inbetween
          bin_len_list <- c(bin_len_list, bin_len) #make list of all bin lengths
        }
        final_bin$bin_length <- unlist(bin_len_list) # add to dataframe
        
        
        bin_len_list <- as.numeric(final_bin$bin_length) #new lists since dataframe has been reordered
        count_list <- as.numeric(final_bin$count)
        
        
        bin_len_cum <- list()
        cum_mut_list <- list()
        cum_mut_list[[1]] <- count_list[[1]]  #start list with first bin length and first mutation count
        bin_len_cum[[1]] <- bin_len_list[[1]]
        
        for (u in 2:length(count_list)) {     #repeat for all other bins
          cum_mut_list[[u]] <- cum_mut_list[[u-1]] + count_list[[u]]  #add all calculated bin lengths and mutation counts to value one row above
          bin_len_cum[[u]] <- bin_len_cum[[u-1]] + bin_len_list[[u]]
        }
        
        final_bin$Cummulative_Bin_Length <- unlist(bin_len_cum)
        final_bin$Cummulative_Mutations <- unlist(cum_mut_list)
        
        #######################################
        ## Apply Cutoff of BP Length if Needed
        #######################################
        
        final_bin <- final_bin[final_bin$Cummulative_Bin_Length <= panel_length,] ##apply bp length cutoff if needed
        
        final_bin <- final_bin[,c(1:4,7:8)]
        
        # add gene names
        gene_name_list <- list()
        
        include_genes <- req(input$radio)
        
        if (include_genes == 1){
          for (j in 1:nrow(final_bin)){
            data_subset <- subset(data, chr == final_bin$chromosome[[j]])
            vec <- final_bin$upperbound[[j]]:final_bin$lowerbound[[j]]
            int <- intersect(data_subset$pos, vec)
            gene_list_temp <- c()
            for (i in int){
              int_sub <- subset(data_subset, pos == i)
              gene_list_temp <- c(gene_list_temp, int_sub$gene)
            }
            name <- paste(shQuote(unique(gene_list_temp)), collapse=", ")
            gene_name_list <- c(gene_name_list, name)
          }
          final_bin$gene <- unlist(gene_name_list)
          colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION.COUNT", "CUMULATIVE.PANEL.LENGTH",
                                   "CUMULATIVE MUTATIONS CAPTURED", "GENE")
          
          myPlot(ggplotly(ggplot(final_bin, aes(x = CUMULATIVE.PANEL.LENGTH, y = MUTATION.COUNT))+
                            geom_point() +
                            geom_line() +
                            ylab("Mutations Captured per Amplicon") +
                            xlab("Sequencing Panel Length") +
                            theme_classic()+
                            theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5, size = 10), 
                                  axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10))))
          
          colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION COUNT", "CUMULATIVE PANEL LENGTH",
                                   "CUMULATIVE MUTATIONS CAPTURED", "GENE")
          final_bin <- final_bin[,c(3,1:2,4,7,5:6)]
          rownames(final_bin) <- 1:nrow(final_bin)
          
          gene_freq <- data.frame(ftable(final_bin$GENE)) #make frequency table for each position
          gene_freq <- gene_freq[order(-gene_freq$Freq), ]
          
          
          list_gene_lengths <- list()
          unique_genes <- gene_freq[,1]
          for (i in 1:length(unique_genes)){
            temp <- paste(unique_genes[[i]], ":", nrow(subset(final_bin, GENE == unique_genes[[i]])))
            list_gene_lengths <- c(list_gene_lengths, temp)
          }
          
          text(paste("Total number of amplicons: ", nrow(final_bin), ", Average number of mutations in amplicon: ", mean(final_bin[,4]),
               "Genes: ", toString(list_gene_lengths)))
        }
        

        if (include_genes == 2){
        colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION.COUNT", "CUMULATIVE.PANEL.LENGTH",
                                 "CUMULATIVE MUTATIONS CAPTURED")
        
        myPlot(ggplotly(ggplot(final_bin, aes(x = CUMULATIVE.PANEL.LENGTH, y = MUTATION.COUNT))+
                            geom_point() +
                            geom_line() +
                            ylab("Mutations Captured per Amplicon") +
                            xlab("Sequencing Panel Length") +
                            theme_classic()+
                            theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5, size = 10), 
                                  axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10))))
        
        colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION COUNT", "CUMULATIVE PANEL LENGTH",
                                 "CUMULATIVE MUTATIONS CAPTURED")
        final_bin <- final_bin[,c(3,1:2,4:6)]
        rownames(final_bin) <- 1:nrow(final_bin)
        
        text(paste("Total number of amplicons: ", nrow(final_bin), ", Average number of mutations in amplicon: ", mean(final_bin[,4])))
        
        }
        
        myData(final_bin)
        myFile(final_bin)


        

        
        
    })
        
        output$myPlot = renderPlotly({
            myPlot()
        })
        output$myData = renderDataTable({
            myData()
        })
        output$text = renderText({
          text()
        })
        
        output$myFile = downloadHandler(filename = function (){paste0("hotSPOT_sequencing_panel", Sys.Date(), ".csv", sep ="")},
          content = function(file){write.csv(myFile(), file)}, contentType = "text/csv")
        
        myPlot2 = reactiveVal()
        myData2 = reactiveVal()
        myFile2 = reactiveVal()
        text2 = reactiveVal()
        
    
    
        observeEvent(input$goButton2, {
        
        data2 <- req(input$data2)
        data <- read.csv(data2$datapath)
        
        panel_length = req(input$panel_length2) 
        amp_length = req(input$amp_length2)
        num_amps = req(input$num_amps)
        
        ######### amplicon based binning algorithm
        
        library(dplyr)
        library(plyr)
        library(hash)
        library(rlist)
        library(R.utils)
        library(ggplot2)
        library(plotly)
        
        
        #######################################################################
        ## run forward binning method first to establish an appropriate cutoff
        #######################################################################
        ######################################################################
        ## dataset must be adjusted to have two columns: "chr" with chromosome positions
        ## and "pos" with bp position ^^^
        ######################################################################
        
        pos <- data$pos
        pos_freq <- data.frame(ftable(pos)) #make frequency table for each position
        
        data <- merge(data, pos_freq, by = "pos") #merge frequency table with original df
        data <- unique(data) #keep only unique positions
        
        
        amplicon_finder <- function(data, chr){
          pos <- data$pos
          
          keys <- make.keys(data$pos)
          dict <- hash(keys = keys, values = data$Freq) #create dictionary for mutation locations and mutation frequency
          
          make_bins <- sort(unique(pos))
          bins <- data.frame()
          while (length(make_bins) > 1){
            temp_bins <- make_bins
            i <- 1
            start <- make_bins[1]
            p_dif <- 1
            go = "yes"
            while (go == "yes"){
              p1 <- temp_bins[i]
              p_next <- temp_bins[i + 1]
              i <- i + 1
              if (i <= length(make_bins)){
                p_dif <- p_next - p1
              }
              if (i > length(make_bins)){
                p_dif = amp_length + 100
                go = "no"
              }
              if (p_dif > amp_length){
                p_next <- p1
                go = "no"
              }
            }
            end <- p_next
            vec <- as.list(start:end)
            keep <- setdiff(make_bins, vec)
            make_bins <- keep
            chr = chr
            row <- data.frame("lowerbound" = start, "upperbound" = end, "chromosome" = chr)
            bins <- rbind(bins, row)
          }
          if (length(make_bins) == 1){
            row <- data.frame("lowerbound" = make_bins, "upperbound" = make_bins, "chromosome" = chr)
            bins <- rbind(bins, row)
          }
          
          mut_locations <- data.frame()
          possible_bins <- data.frame()
          
          for (i in 1:nrow(bins)){
            up_mut <- bins$upperbound[[i]]
            low_mut <- bins$lowerbound[[i]]
            vec <- up_mut:low_mut
            dif <- length(vec) - amp_length
            if (dif <= 0){ #############FIX HERE#######################
              up <- bins$upperbound[[i]]
              low <- bins$lowerbound[[i]]
              vec <- up:low
              mutations <- intersect(vec, pos) #find all mutations in total bin region
              mutations_keys <- make.keys(mutations)
              count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
              mut_all <- c()
              for (m in 1:length(mutations)){
                mut_all <- c(mut_all, rep(mutations[[m]], as.vector(values(dict, keys = mutations_keys))[[m]]))
              }
              weighted_mid <- round(mean(mut_all))
              weighted_mid_max <- (ceiling(mean(vec)) + (amp_length/2))
              weighted_mid_min <- (floor(mean(vec)) - (amp_length/2))
              if (weighted_mid < weighted_mid_min){final_mid <- weighted_mid_min}
              if (weighted_mid > weighted_mid_max){final_mid <- weighted_mid_max}
              if (weighted_mid > weighted_mid_min & weighted_mid < weighted_mid_max){final_mid <- weighted_mid}
              up <- final_mid + ceiling((amp_length - 1) / 2)
              low <- final_mid - floor((amp_length - 1) / 2)
              vec <- up:low
              mutations <- intersect(vec, pos) #find all mutations in total bin region
              mutations_keys <- make.keys(mutations)
              count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
              row <- data.frame("lowerbound" = low, "upperbound" = up, 
                                "chromosome" = bins$chromosome[[i]], "count" = count,"id" = "x")
              possible_bins <- rbind(possible_bins, row)
            }
            else{
              up_mut <- bins$upperbound[[i]]
              low_mut <- bins$lowerbound[[i]]
              vec <- up_mut:low_mut
              mut_list <- intersect(vec, pos)
              id = paste(as.character(bins$chromosome[[i]]), as.character(i), sep = "-")
              b_list <- list()
              for (mut in mut_list){
                bin1 <- (mut - amp_length):(mut - 1)
                b_list[[1]] <- bin1
                bin2 <- (mut - (amp_length - 1)):(mut)
                b_list[[2]] <- bin2
                bin3 <- (mut + amp_length):(mut + 1)
                b_list[[3]] <- bin3
                bin4 <- (mut + (amp_length - 1)):(mut)
                b_list[[4]] <- bin4
                bin5 <- (mut - ceiling(amp_length / 2)):(mut - floor(amp_length / 2))
                b_list[[5]] <- bin5
                
                for (b in b_list){
                  low <- min(b)
                  up <- max(b)
                  mutations <- intersect(b, pos)
                  if (length(mutations > 0)){
                    mutations_keys <- make.keys(mutations)
                    count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                    row <- data.frame("lowerbound" = low, "upperbound" = up, 
                                      "chromosome" = bins$chromosome[[i]], "count" = count,"id" = id)
                    possible_bins <- rbind(possible_bins, row)
                  }
                }
              }
              
            }
          }
          return(possible_bins)
        }
        
        ##################################################################################
        
        
        chrom_list <- unique(as.list(data$chr)) # list of all chromosomes found in this dataset
        all_pos_bins <- data.frame()
        
        for (chrom in chrom_list) {
          chromosome <- subset(data, chr == chrom)
          chromsort <- chromosome[order(chromosome$pos), ] #order df from lowest to highest position
          chrom_bins <- amplicon_finder(chromsort, chrom) #run through above function
          all_pos_bins <- rbind(all_pos_bins, chrom_bins) #add all bins together from all chromosomes
        }
        
        ###################################################################################
        ordered_bin <- all_pos_bins[order(-all_pos_bins$count), ]
        
        big_spots <- ordered_bin[!ordered_bin$id == "x", ]
        unique_spots <- unique(as.list(big_spots$id))
        
        keys <- make.keys(data$pos)
        dict <- hash(keys = keys, values = data$Freq)
        possible_bins <- data.frame()
        for (spot in unique_spots){
          spot_sub <- subset(big_spots, id == spot)
          big_bin_all <- list()
          while (nrow(spot_sub) > 0){
            row <- data.frame()
            new_sub <- data.frame()
            big_spot <- spot_sub[1,]
            big_bin <- big_spot$upperbound:big_spot$lowerbound
            big_bin_all <- c(big_bin_all, big_bin)
            possible_bins <- rbind(possible_bins, big_spot)
            if (nrow(spot_sub) > 1){
              for (k in 2:nrow(spot_sub)){
                row <- data.frame()
                test_bin <- spot_sub$upperbound[k]:spot_sub$lowerbound[k]
                overlap <- intersect(big_bin_all, test_bin)
                if (length(overlap > 0)){
                  unique_bin <- setdiff(test_bin, overlap)
                  mutations <- intersect(unique_bin, pos) #find all mutations in total bin region
                  if (length(mutations) > 0){
                    mutations_keys <- make.keys(mutations)
                    new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                    row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin), 
                                      "chromosome" = spot_sub$chromosome[1], "count" = new_count,"id" = spot_sub$id[1])
                  }
                }
                if (length(overlap) <= 0) {row <- spot_sub[k,]}
                new_sub <- rbind(new_sub, row)
              }
            }
            spot_sub <- new_sub
            if (nrow(spot_sub > 0)){spot_sub <- spot_sub[order(-spot_sub$count), ]}
          }
        }
        
        
        final_bin <- rbind(subset(ordered_bin, id == "x"), possible_bins)
        final_bin <- final_bin[order(-final_bin$count), ]
        
        
        
        ##############################################################
        
        
        bin_len_list <- list()
        for (y in 1:nrow(final_bin)) {      #calculating bin length using upper and lower bound positions
          upper = final_bin$upperbound[[y]]
          lower = final_bin$lowerbound[[y]]
          bin_len <- upper - lower + 1# add 1 to include both upper and lowerbound plus all bp inbetween
          bin_len_list <- c(bin_len_list, bin_len) #make list of all bin lengths
        }
        final_bin$bin_length <- unlist(bin_len_list) # add to dataframe
        
        
        bin_len_list <- as.numeric(final_bin$bin_length) #new lists since dataframe has been reordered
        count_list <- as.numeric(final_bin$count)
        
        
        bin_len_cum <- list()
        cum_mut_list <- list()
        cum_mut_list[[1]] <- count_list[[1]]  #start list with first bin length and first mutation count
        bin_len_cum[[1]] <- bin_len_list[[1]]
        
        for (u in 2:length(count_list)) {     #repeat for all other bins
          cum_mut_list[[u]] <- cum_mut_list[[u-1]] + count_list[[u]]  #add all calculated bin lengths and mutation counts to value one row above
          bin_len_cum[[u]] <- bin_len_cum[[u-1]] + bin_len_list[[u]]
        }
        
        final_bin$Cummulative_Bin_Length <- unlist(bin_len_cum)
        final_bin$Cummulative_Mutations <- unlist(cum_mut_list)
        
        # find cutoff value
        length <- nrow(final_bin[final_bin$Cummulative_Bin_Length <= panel_length,])
        cutoff <- final_bin$count[[length]]
        
        ###################################################################################
        rm(list= ls()[!(ls() %in% c("data", "cutoff", "num_amps", "panel_length", "pos", "amp_length"))])
        
        
        amplicon_finder <- function(data, chr){
          pos <- data$pos
          
          keys <- make.keys(data$pos)
          dict <- hash(keys = keys, values = data$Freq) #create dictionary for mutation locations and mutation frequency
          
          make_bins <- sort(unique(pos))
          bins <- data.frame()
          while (length(make_bins) > 1){
            temp_bins <- make_bins
            i <- 1
            start <- make_bins[1]
            p_dif <- 1
            go = "yes"
            while (go == "yes"){
              p1 <- temp_bins[i]
              p_next <- temp_bins[i + 1]
              i <- i + 1
              if (i <= length(make_bins)){
                p_dif <- p_next - p1
              }
              if (i > length(make_bins)){
                p_dif = amp_length + 100
                go = "no"
              }
              if (p_dif > amp_length){
                p_next <- p1
                go = "no"
              }
            }
            end <- p_next
            vec <- as.list(start:end)
            keep <- setdiff(make_bins, vec)
            make_bins <- keep
            chr = chr
            row <- data.frame("lowerbound" = start, "upperbound" = end, "chromosome" = chr)
            bins <- rbind(bins, row)
          }
          if (length(make_bins) == 1){
            row <- data.frame("lowerbound" = make_bins, "upperbound" = make_bins, "chromosome" = chr)
            bins <- rbind(bins, row)
          }
          
          mut_locations <- data.frame()
          possible_bins <- data.frame()
          
          for (i in 1:nrow(bins)){
            up_mut <- bins$upperbound[[i]]
            low_mut <- bins$lowerbound[[i]]
            vec <- up_mut:low_mut
            dif <- length(vec) - amp_length
            if (dif <= 0){ #############FIX HERE#######################
              up <- bins$upperbound[[i]]
              low <- bins$lowerbound[[i]]
              vec <- up:low
              mutations <- intersect(vec, pos) #find all mutations in total bin region
              mutations_keys <- make.keys(mutations)
              count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
              mut_all <- c()
              for (m in 1:length(mutations)){
                mut_all <- c(mut_all, rep(mutations[[m]], as.vector(values(dict, keys = mutations_keys))[[m]]))
              }
              weighted_mid <- round(mean(mut_all))
              weighted_mid_max <- (ceiling(mean(vec)) + (amp_length/2))
              weighted_mid_min <- (floor(mean(vec)) - (amp_length/2))
              if (weighted_mid < weighted_mid_min){final_mid <- weighted_mid_min}
              if (weighted_mid > weighted_mid_max){final_mid <- weighted_mid_max}
              if (weighted_mid > weighted_mid_min & weighted_mid < weighted_mid_max){final_mid <- weighted_mid}
              up <- final_mid + ceiling((amp_length - 1) / 2)
              low <- final_mid - floor((amp_length - 1) / 2)
              vec <- up:low
              mutations <- intersect(vec, pos) #find all mutations in total bin region
              mutations_keys <- make.keys(mutations)
              count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
              row <- data.frame("lowerbound" = low, "upperbound" = up, 
                                "chromosome" = bins$chromosome[[i]], "count" = count,"id" = "x", "mut_lowerbound" = low, "mut_upperbound" = up)
              possible_bins <- rbind(possible_bins, row)
            }
            else{
              up_mut <- bins$upperbound[[i]]
              low_mut <- bins$lowerbound[[i]]
              vec <- up_mut:low_mut
              mut_list <- intersect(vec, pos)
              id = paste(as.character(bins$chromosome[[i]]), as.character(i), sep = "-")
              b_list <- list()
              for (mut in mut_list){
                bin1 <- (mut - amp_length):(mut - 1)
                b_list[[1]] <- bin1
                bin2 <- (mut - (amp_length - 1)):(mut)
                b_list[[2]] <- bin2
                bin3 <- (mut + amp_length):(mut + 1)
                b_list[[3]] <- bin3
                bin4 <- (mut + (amp_length - 1)):(mut)
                b_list[[4]] <- bin4
                bin5 <- (mut - ceiling(amp_length / 2)):(mut - floor(amp_length / 2))
                b_list[[5]] <- bin5
                
                for (b in b_list){
                  low <- min(b)
                  up <- max(b)
                  mutations <- intersect(b, pos)
                  if (length(mutations > 0)){
                    mutations_keys <- make.keys(mutations)
                    count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                    row <- data.frame("lowerbound" = low, "upperbound" = up, 
                                      "chromosome" = bins$chromosome[[i]], "count" = count,"id" = id, "mut_lowerbound" = low_mut, "mut_upperbound" = up_mut)
                    possible_bins <- rbind(possible_bins, row)
                  }
                }
              }
              
            }
          }
          return(possible_bins)
        }
        
        
        
        ##################################################################################
        
        
        chrom_list <- unique(as.list(data$chr)) # list of all chromosomes found in this dataset
        all_pos_bins <- data.frame()
        
        for (chrom in chrom_list) {
          chromosome <- subset(data, chr == chrom)
          chromsort <- chromosome[order(chromosome$pos), ] #order df from lowest to highest position
          chrom_bins <- amplicon_finder(chromsort, chrom) #run through above function
          all_pos_bins <- rbind(all_pos_bins, chrom_bins) #add all bins together from all chromosomes
        }
        ordered_bin <- all_pos_bins[order(-all_pos_bins$count), ]
        
        big_spots <- ordered_bin[!ordered_bin$id == "x", ]
        unique_spots <- unique(as.list(big_spots$id))
        
        keys <- make.keys(data$pos)
        dict <- hash(keys = keys, values = data$Freq)
        
        
        possible_bins <- data.frame()
        maybe_possible_bins <- data.frame()
        y = 1
        new_pos <- pos
        
        ####################################################################################
        
        new_big_spots <- data.frame()
        maybe_new_big_spots <- data.frame()
        too_big_split <- data.frame()
        new_big_split <- data.frame()
        
        for (spot in unique_spots){ 
          all_vec <- list()
          all_bin <- list()
          spot_sub <- subset(big_spots, id == spot) #subset df by spot id
          num_amplicons <- ceiling(length(spot_sub$mut_lowerbound[[1]]:spot_sub$mut_upperbound[[1]]) / amp_length)
          if (num_amplicons <= num_amps){
            new_big_spots <- rbind(new_big_spots, spot_sub)
          }
          if (num_amplicons > num_amps){
            spot_sub <- subset(spot_sub, count >= cutoff)
            if (nrow(spot_sub) > 0){
              for (i in 1:nrow(spot_sub)){
                vec <- spot_sub$upperbound[[i]]:spot_sub$lowerbound[[i]]
                all_vec <- c(all_vec, vec)
                all_bin[[i]] <- vec
              }
              new_split <- as.data.frame(seqToIntervals(unique(all_vec)))
              too_big_split <- new_split
              
              for (j in 1:nrow(new_split)){
                split_vec <- new_split$to[[j]]:new_split$from[[j]]
                all_mutations <- intersect(split_vec, pos)
                new_id <- paste(as.character(spot), as.character(j), sep = "-")
                for (v in all_bin){
                  split_overlap <- intersect(v, split_vec)
                  if (length(split_overlap) == amp_length){
                    mutations <- intersect(v, pos) #find all mutations in total bin region
                    mutations_keys <- make.keys(mutations)
                    new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                    row <- data.frame("lowerbound" = min(v), "upperbound" = max(v), "count" = new_count,
                                      "chromosome" = spot_sub$chromosome[[1]], "id" = new_id, "mut_lowerbound" = min(all_mutations), 
                                      "mut_upperbound" = max(all_mutations))
                    maybe_new_big_spots <- rbind(maybe_new_big_spots, row)
                  }
                }
              }
              big_amp_ids <- unique(maybe_new_big_spots$id)
              
              for (p in big_amp_ids){
                amp_list <- list()
                
                amp_spot_sub <- subset(maybe_new_big_spots, id == p)
                amp_test <- ceiling(length(amp_spot_sub$mut_lowerbound[[1]]:amp_spot_sub$mut_upperbound[[1]]) / 125)
                if (amp_test <= num_amps){new_big_spots <- rbind(new_big_spots, amp_spot_sub)}
                if (amp_test > num_amps){
                  too_big_split <- data.frame("to" = min(amp_spot_sub$mut_lowerbound), "from" = max(amp_spot_sub$mut_upperbound))
                  while (nrow(too_big_split) > 0){
                    new_amp_spot_sub <- data.frame()
                    too_big_list <- c()
                    for (i in 1:nrow(too_big_split)){
                      too_big_list <- c(too_big_list, too_big_split$from[[i]]:too_big_split$to[[i]])
                    }
                    for (v in 1:nrow(amp_spot_sub)){
                      vec <- amp_spot_sub$lowerbound[[v]]:amp_spot_sub$upperbound[[v]]
                      ov <- intersect(vec, too_big_list)
                      if (length(ov) == amp_length){new_amp_spot_sub <- rbind(new_amp_spot_sub, amp_spot_sub[v,])}
                    }
                    new_amp_spot_sub$mut_upperbound <- max(too_big_list)
                    new_amp_spot_sub$mut_lowerbound <- min(too_big_list)
                    all_big_vec <- list()
                    for (m in 1:nrow(new_amp_spot_sub)){
                      vec <- new_amp_spot_sub$upperbound[[m]]:new_amp_spot_sub$lowerbound[[m]]
                      all_big_vec <- c(all_big_vec, vec)
                      vec_freq <- data.frame(ftable(unlist(all_big_vec)))
                    }
                    low_point <- list()
                    for (x in 1:nrow(vec_freq)){
                      freq_list <- vec_freq$Freq
                      pos_list <- as.numeric(levels(factor(vec_freq$Var1)))
                      f <- freq_list[[x]]
                      minimum <- min(freq_list)
                      if (f == minimum){
                        low_point <- c(low_point, pos_list[[x]])
                      }
                    }
                    low_point_split <- as.data.frame(seqToIntervals(low_point))
                    list_length <- list()
                    for (low in 1:nrow(low_point_split)){
                      length <- length(low_point_split$to[[low]]:low_point_split$from[[low]])
                      list_length <- c(list_length, length)
                    }
                    long_point <- which.max(list_length) 
                    break_point <- ceiling(median(low_point_split$to[[long_point]]:low_point_split$from[[long_point]]))
                    all_big_vec <- vec_freq$Var1
                    all_big_vec <- as.numeric(levels(factor(all_big_vec)))
                    all_big_vec <- all_big_vec[all_big_vec != break_point]
                    test_big_split <- as.data.frame(seqToIntervals(unique(all_big_vec)))
                    too_big_split <- data.frame()
                    for (u in 1:nrow(test_big_split)){
                      vec <- test_big_split$to[[u]]:test_big_split$from[[u]]
                      len <- length(vec)
                      if (len <= (num_amps * amp_length)){
                        new_big_split <- rbind(new_big_split, test_big_split[u,])
                      }
                      if (len > (num_amps * amp_length)){
                        too_big_split <- rbind(too_big_split, test_big_split[u,])
                      }
                    }
                  }
                  for (j in 1:nrow(new_big_split)){
                    split_big_vec <- new_big_split$to[[j]]:new_big_split$from[[j]]
                    all_mutations <- intersect(split_big_vec, pos)
                    new_id <- paste(as.character(new_id), as.character(j), sep = "-")
                    for (v in all_bin){
                      split_overlap <- intersect(v, split_big_vec)
                      if (length(split_overlap) >= amp_length){
                        mutations <- intersect(v, pos) #find all mutations in total bin region
                        mutations_keys <- make.keys(mutations)
                        new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                        row <- data.frame("lowerbound" = min(v), "upperbound" = max(v), "count" = new_count,
                                          "chromosome" = spot_sub$chromosome[[1]], "id" = new_id, "mut_lowerbound" = min(all_mutations), 
                                          "mut_upperbound" = max(all_mutations))
                        new_big_spots <- rbind(new_big_spots, row)
                      }
                    }
                  }
                }
              }
            }
          }
        }
        ########################################################################################################
        
        all_best_combo <- list()
        fin_possible_bins <- data.frame()
        
        new_unique_spots <- unique(as.list(new_big_spots$id))
        for (spot in new_unique_spots){ 
          spot_sub <- subset(new_big_spots, id == spot) #subset df by spot id
          spot_sub <- subset(spot_sub, count >= cutoff)
          spot_sub <- unique(spot_sub)
          if (nrow(spot_sub) == 1){
            spot_sub <- spot_sub[,1:5]
            fin_possible_bins <- rbind(fin_possible_bins, spot_sub)
          }
          if (nrow(spot_sub) > 1){
            vec_list <- list()
            for (k in 1:nrow(spot_sub)){
              upper <- spot_sub$upperbound[[k]]
              lower <- spot_sub$lowerbound[[k]]
              vec <- upper:lower
              vec_list[[k]] <- vec
            }
            num_amplicons <- ceiling(length(spot_sub$mut_lowerbound[[1]]:spot_sub$mut_upperbound[[1]]) / amp_length)
            vec_list_num <- 1:length(vec_list)
            if (num_amplicons == 1){
              pos_comb <- expand.grid(vec_list_num)
            }
            if (num_amplicons == 2){
              pos_comb <- expand.grid(vec_list_num, vec_list_num)
            }
            if (num_amplicons == 3){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 4){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 5){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 6){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 7){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 8){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 9){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            if (num_amplicons == 10){
              pos_comb <- expand.grid(vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num, vec_list_num)
            }
            all_count <- list()
            fin_sum <- list()
            int_sum <- list()
            list_pos_bins <- list()
            for (c in 1:nrow(pos_comb)){
              com <- pos_comb[c,]
              un_com <- unique(unlist(com))
              if (length(un_com) < length(com)){
                pos_comb <- pos_comb[-c,]
              }
            }
            for (j in 1:nrow(pos_comb)){
              vec_comb <- list()
              temp_df <- data.frame()
              for (i in 1:ncol(pos_comb)){
                vec <- unlist(vec_list[pos_comb[j,i]])
                mutations <- intersect(vec, new_pos) #find all mutations in total bin region
                mutations_keys <- make.keys(mutations)
                count <- sum(as.vector(values(dict, keys = mutations_keys)))
                row <- data.frame("lowerbound" = min(vec), "upperbound" = max(vec),"count" = count)
                temp_df <- rbind(temp_df, row)
              }
              maybe_bins_1 <- temp_df[order(-temp_df$count), ]
              maybe_bins_2 <- data.frame(maybe_bins_1) #
              maybe_bins_2$chromosome <- spot_sub$chromosome[[1]]
              maybe_bins_2$id <- spot_sub$id[[1]]
              list_pos_bins[[j]] <- maybe_bins_2
              for (x in 1:nrow(maybe_bins_2)){
                sum_count <- list()
                maybe_bins <- maybe_bins_2
                maybe_bins[1,] <- maybe_bins_1[x,]
                maybe_bins[x,] <- maybe_bins_1[1,]
                big_bin_all <- list()
                possible_bins <- data.frame()
                while (nrow(maybe_bins) > 0){
                  row <- data.frame()
                  new_bin <- data.frame()
                  big_spot <- maybe_bins[1,]
                  big_bin <- big_spot$upperbound:big_spot$lowerbound
                  big_bin_all <- c(big_bin_all, big_bin)
                  possible_bins <- rbind(possible_bins, big_spot)
                  if (nrow(maybe_bins) > 1){
                    for (k in 2:nrow(maybe_bins)){
                      row <- data.frame()
                      test_bin <- maybe_bins$upperbound[k]:maybe_bins$lowerbound[k]
                      overlap <- intersect(big_bin_all, test_bin)
                      if (length(overlap) > 0){
                        unique_bin <- setdiff(test_bin, overlap)
                        mutations <- intersect(unique_bin, new_pos) #find all mutations in total bin region
                        if (length(mutations) > 0){
                          mutations_keys <- make.keys(mutations)
                          new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                          row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin), "count" = new_count,
                                            "chromosome" = spot_sub$chromosome[[1]], "id" = spot_sub$id[[1]])
                        }
                      }
                      if (length(overlap) <= 0) {row <- maybe_bins[k,]}
                      new_bin <- rbind(new_bin, row)
                    }
                  }
                  maybe_bins <- new_bin
                  if (nrow(maybe_bins > 0)){maybe_bins <- maybe_bins[order(-maybe_bins$count), ]}
                }
                testing_bins <- subset(possible_bins, count >= cutoff)
                sum <- sum(unlist(testing_bins$count))
                sum_count <- c(sum_count, sum)
              }
              fin_sum <- c(fin_sum, max(unlist(sum_count)))
              int_sum <- c(int_sum, which.max(sum_count))
            }
            best_combo <- which.max(fin_sum)
            best_one <- as.numeric(best_combo)
            bins <- list_pos_bins[[best_one]]
            best_int <- int_sum[[best_one]]
            maybe_bins <- bins
            maybe_bins[1,] <- bins[as.numeric(best_int),]
            maybe_bins[as.numeric(best_int),] <- bins[1,]
            big_bin_all <- list()
            while (nrow(maybe_bins) > 0){
              row <- data.frame()
              new_bin <- data.frame()
              big_spot <- maybe_bins[1,]
              big_bin <- big_spot$upperbound:big_spot$lowerbound
              mutations <- intersect(big_bin, new_pos)
              mutations_keys <- make.keys(mutations)
              new_count <- sum(as.vector(values(dict, keys = mutations_keys))) 
              big_bin_all <- c(big_bin_all, big_bin)
              fin_possible_bins <- rbind(fin_possible_bins, big_spot)
              if (nrow(maybe_bins) > 1){
                for (k in 2:nrow(maybe_bins)){
                  row <- data.frame()
                  test_bin <- maybe_bins$upperbound[k]:maybe_bins$lowerbound[k]
                  overlap <- intersect(big_bin_all, test_bin)
                  if (length(overlap) > 0){
                    unique_bin <- setdiff(test_bin, overlap)
                    mutations <- intersect(unique_bin, new_pos) #find all mutations in total bin region
                    if (length(mutations) > 0){
                      mutations_keys <- make.keys(mutations)
                      new_count <- sum(as.vector(values(dict, keys = mutations_keys))) #get total mutations for total bin region
                      row <- data.frame("lowerbound" = min(test_bin), "upperbound" = max(test_bin), "count" = new_count, 
                                        "chromosome" = spot_sub$chromosome[[1]], id = spot_sub$id[[1]])
                    }
                  }
                  if (length(overlap) <= 0) {row <- maybe_bins[k,]}
                  new_bin <- rbind(new_bin, row)
                }
              }
              maybe_bins <- new_bin
              if (nrow(maybe_bins > 0)){maybe_bins <- maybe_bins[order(-maybe_bins$count), ]}
            }
            
          }
        }
        
        
        best_bin_all <- fin_possible_bins
        ordered_bin <- ordered_bin[-c(6:7)]
        
        final_bin <- rbind(subset(ordered_bin, id == "x"), best_bin_all)
        final_bin <- final_bin[order(-final_bin$count), ]
        
        ##############################################################
        # reorder to favor amplicons with mutations in center
        
        cutoff_point <- final_bin$count[[round(panel_length / amp_length)]]
        
        subset_final_bin <- subset(final_bin, count == cutoff_point)
        dif_list <- list()
        for (i in nrow(subset_final_bin)){
          whole_region <- subset_final_bin$lowerbound[[i]]:subset_final_bin$upperbound[[i]]
          bin_mutations <- intersect(pos, whole_region)
          mid_diff <- abs((max(whole_region) - max(bin_mutations)) - (min(bin_mutations) - min(whole_region)))
          dif_list <- c(dif_list, mid_diff)
        }
        subset_final_bin$dif <- unlist(dif_list)
        subset_final_bin <- subset_final_bin[order(-subset_final_bin$dif), ]
        subset_final_bin <- subset_final_bin[,-ncol(subset_final_bin)]
        
        df_range_min <- min(which(final_bin$count == cutoff_point))
        df_range_max <- max(which(final_bin$count == cutoff_point))
        
        final_bin[df_range_min:df_range_max,] <- subset_final_bin[1:nrow(subset_final_bin),]
        
        
        bin_len_list <- list()
        for (y in 1:nrow(final_bin)) {      #calculating bin length using upper and lower bound positions
          upper = final_bin$upperbound[[y]]
          lower = final_bin$lowerbound[[y]]
          bin_len <- upper - lower + 1# add 1 to include both upper and lowerbound plus all bp inbetween
          bin_len_list <- c(bin_len_list, bin_len) #make list of all bin lengths
        }
        final_bin$bin_length <- unlist(bin_len_list) # add to dataframe
        
        
        bin_len_list <- as.numeric(final_bin$bin_length) #new lists since dataframe has been reordered
        count_list <- as.numeric(final_bin$count)
        
        
        bin_len_cum <- list()
        cum_mut_list <- list()
        cum_mut_list[[1]] <- count_list[[1]]  #start list with first bin length and first mutation count
        bin_len_cum[[1]] <- bin_len_list[[1]]
        
        for (u in 2:length(count_list)) {     #repeat for all other bins
          cum_mut_list[[u]] <- cum_mut_list[[u-1]] + count_list[[u]]  #add all calculated bin lengths and mutation counts to value one row above
          bin_len_cum[[u]] <- bin_len_cum[[u-1]] + bin_len_list[[u]]
        }
        
        final_bin$Cummulative_Bin_Length <- unlist(bin_len_cum)
        final_bin$Cummulative_Mutations <- unlist(cum_mut_list)
        
        #######################################
        ## Apply Cutoff of BP Length if Needed
        #######################################
        #bins_final_new <- unique(bins_ordered_new)
        final_bin <- final_bin[final_bin$Cummulative_Bin_Length <= panel_length,] ##apply bp length cutoff if needed
        final_bin <- final_bin[,c(1:4,7:8)]
        # add gene names
        gene_name_list <- list()
        
        include_genes <- req(input$radio)
        
        if (include_genes == 1){
          for (j in 1:nrow(final_bin)){
            data_subset <- subset(data, chr == final_bin$chromosome[[j]])
            vec <- final_bin$upperbound[[j]]:final_bin$lowerbound[[j]]
            int <- intersect(data_subset$pos, vec)
            gene_list_temp <- c()
            for (i in int){
              int_sub <- subset(data_subset, pos == i)
              gene_list_temp <- c(gene_list_temp, int_sub$gene)
            }
            name <- paste(shQuote(unique(gene_list_temp)), collapse=", ")
            gene_name_list <- c(gene_name_list, name)
          }
          final_bin$gene <- unlist(gene_name_list)
          colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION.COUNT", "CUMULATIVE.PANEL.LENGTH",
                                   "CUMULATIVE MUTATIONS CAPTURED", "GENE")
          
          myPlot2(ggplotly(ggplot(final_bin, aes(x = CUMULATIVE.PANEL.LENGTH, y = MUTATION.COUNT))+
                            geom_point() +
                            geom_line() +
                            ylab("Mutations Captured per Amplicon") +
                            xlab("Sequencing Panel Length") +
                            theme_classic()+
                            theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5, size = 10), 
                                  axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10))))
          
          colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION COUNT", "CUMULATIVE PANEL LENGTH",
                                   "CUMULATIVE MUTATIONS CAPTURED", "GENE")
          final_bin <- final_bin[,c(3,1:2,4,7,5:6)]
          rownames(final_bin) <- 1:nrow(final_bin)
        }
        
        
        if (include_genes == 2){
          colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION.COUNT", "CUMULATIVE.PANEL.LENGTH",
                                   "CUMULATIVE MUTATIONS CAPTURED")
          
          myPlot(ggplotly(ggplot(final_bin, aes(x = CUMULATIVE.PANEL.LENGTH, y = MUTATION.COUNT))+
                            geom_point() +
                            geom_line() +
                            ylab("Mutations Captured per Amplicon") +
                            xlab("Sequencing Panel Length") +
                            theme_classic()+
                            theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=0.5, size = 10), 
                                  axis.title.y = element_text(size = 10), axis.text.y = element_text(size = 10), title = element_text(size = 10))))
          
          colnames(final_bin) <- c("LOWERBOUND", "UPPERBOUND", "CHROMOSOME", "MUTATION COUNT", "CUMULATIVE PANEL LENGTH",
                                   "CUMULATIVE MUTATIONS CAPTURED")
          final_bin <- final_bin[,c(3,1:2,4:6)]
          rownames(final_bin) <- 1:nrow(final_bin)
          
        }
        
        myData2(final_bin)
        myFile2(final_bin)
        text2(paste("Total number of amplicons: ", nrow(final_bin), ", Average number of mutations in amplicon: ", mean(final_bin[,4])))
        
        
    })
        output$myPlot2 = renderPlotly({
            myPlot2()
        })
        output$myData2 = renderDataTable({
            myData2()
        })
        output$text2 = renderText({
          text2()
        })
        
        output$myFile2 = downloadHandler(filename = function (){paste0("hotSPOT_sequencing_panel", Sys.Date(), ".csv", sep ="")},
                                        content = function(file){write.csv(myFile2(), file)}, contentType = "text/csv")
}

shinyApp(ui = ui, server = server)
