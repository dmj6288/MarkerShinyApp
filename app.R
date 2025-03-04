

source("baseline.R")

options(shiny.fullstacktrace = TRUE)


# Load Shiny library
library(shiny)

# Define the UI (User Interface)
ui <- fluidPage(
  titlePanel("Yi lab marker generator"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("pvalue", "Maximum Adjusted p-value", min = 0, max = 0.25, value = 0.05, step = 0.01),
      sliderInput("spec", "Minimum Specificity value", min = 0, max = 1, value = 0.65, step = 0.01),
      sliderInput("raw_conc", "Minimum raw transcript concentration", min = 0, max = 1000, value = 100, step = 50),
      sliderInput("scale_conc", "Minimum scaled transcript concentration", min = 0, max = 5, value = 2, step = 0.1),
      sliderInput("FC", "Minimum Fold Change", min = 0, max = 4, value = 1, step = 0.1),
      sliderInput("plotHeight", "Heatmap Height:", min = 600, max = 2400, value = 1200)
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", 
                 
                 div(
                   class = "plot-container",
                   plotOutput("heatmapPlot",  width = "100%", height = "2000px")
                   )),
         tabPanel("Marker lists",
                  div(
                    class = "table-container",
                    DTOutput("gene_table")
                  )
         )
      )
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Update column choices dynamically based on uploaded data
  
  # dataset <- reactive({
  #   req(input$file)  # Ensure a file is uploaded
  #   read.csv(input$file$datapath)
  # })
  # 
  # observeEvent(dataset(), {
  #   updateSelectInput(session, "colSelect", choices = names(dataset()))
  # })
  
  heatmap_294_SW30 <- reactive({
    readRDS("source/sample_wise30_MicroVIS_filtered_pseudobulk_count_matrix_RAW_RN004.rds")
  })
  
  studies_thresholded <- reactive({
    readRDS("source/studies_to_keep_threshold.rds")
  })
  
  DE_tables_with_specificity <- reactive({
    readRDS("source/final_DE_table_with_specificity_RN004.rds")
  })
  
    # Define a reactive expression that returns multiple input values as a list
  user_inputs <- reactive({
  
          specificity_vector = input$spec
          foldchange_vector  = input$FC
          padj_vector        = input$pvalue
          matrix_min         = input$raw_conc
          scale_min          = input$scale_conc
            
          # heatmap_294_SW30           <- heatmap_294_SW30()
          # studies_thresholded        <- studies_thresholded()
          # DE_tables_with_specificity <- DE_tables_with_specificity()
          
          gene_name_lists            <- list()
          gene_name_lists_by_cell    <- list()
          
          parameter_sweep_results    <- c()
          
          for (specificity in specificity_vector){
            
            for (fold_change in foldchange_vector){
              
              for (padj in padj_vector){
                
                for (cell_counts in c("250")){
                  
                  gene_counts_by_cell_type <- c()
                  gene_names_by_cell_type  <- c()
                  
                  for (cell_type in cell_names){
                    
                    region_wise_DE_table <- list()
                    
                    for (region in studies_thresholded()[[cell_counts]][[cell_type]]){
                      
                      current_table <- DE_tables_with_specificity()[[region]][[cell_type]]
                      
                      if (is.null(current_table) == FALSE){
                        
                        filtered_table <- current_table[current_table$Specificity_Index > specificity &
                                                          current_table$log2FoldChange    > fold_change &
                                                       #   current_table$baseMean          > as.numeric(cut_offs[[5]]) &
                                                          current_table$padj              < padj, ]
                        
                        if (length(rownames(filtered_table)) > 0){
                          
                          region_wise_DE_table[[region]] <- rownames(filtered_table) 
                          
                        }
                        
                      }
                      
                      all_genes_for_current_criterion <- Reduce(intersect, region_wise_DE_table)
                      
                    }
                    
                    gene_counts_by_cell_type <- c(gene_counts_by_cell_type, 
                                                  length(all_genes_for_current_criterion), 
                                                  length(studies_thresholded()[[cell_counts]][[cell_type]]))
                    
                    gene_names_by_cell_type  <- c(gene_names_by_cell_type, all_genes_for_current_criterion)
                    
                    gene_name_lists_by_cell[[cell_type]] <- all_genes_for_current_criterion
                    
                  }
                  
                  
                  gene_name_lists[[paste(specificity, fold_change, padj, sep = "_")]] <- gene_names_by_cell_type
                  
                }
                
              }
              
            }
            
          }
          
          list(
            specificity_vector = specificity_vector,
            foldchange_vector  = foldchange_vector,
            padj_vector        = padj_vector,
            matrix_min         = matrix_min,
            scale_min          = scale_min,
            gene_name_lists_by_cell = gene_name_lists_by_cell,
            gene_name_lists = gene_name_lists
          )
  
})
  
  
  marker_clash_reactive <- reactive({
    marker_clash <- c()
    gene_name_lists_by_cell <- user_inputs()$gene_name_lists_by_cell  # ✅ Access reactive value
    
    for (cell1 in cell_names) {
      for (cell2 in cell_names) {
        # print(cell1)
        # print(cell2)
        # print(intersect(gene_name_lists_by_cell[[cell1]], gene_name_lists_by_cell[[cell2]]))
        
        if (cell1 != cell2) {
          marker_clash <- c(marker_clash, intersect(gene_name_lists_by_cell[[cell1]], gene_name_lists_by_cell[[cell2]]))
        }
      }
    }
    
    return(marker_clash)  # ✅ Return the result reactively
  })
  
  
    
  heatmap_data <- reactive({
    
    marker_lists            <- user_inputs()$gene_name_lists  # Extract marker list reactively
    #print(marker_lists)
    gene_name_lists_by_cell <- user_inputs()$gene_name_lists_by_cell
    # specificity_vector      <- user_inputs()$specificity_vector  # ✅ Access reactively
    # foldchange_vector       <- user_inputs()$foldchange_vector
    # padj_vector             <- user_inputs()$padj_vector
    
    all_markers <- unique(unlist(marker_lists))
    
    subset_count_matrix <- heatmap_294_SW30()[rownames(heatmap_294_SW30()) %in% all_markers, ]
    filtered_mat <- subset_count_matrix[apply(subset_count_matrix, 1, mean) >= user_inputs()$matrix_min, ]
    
    if (nrow(filtered_mat) == 0) {
      return(NULL)  # Prevent errors if the matrix is empty
    }
    
    scaled_global_count_matrix <- scale(filtered_mat)
    scaled_global_count_matrix <- scaled_global_count_matrix[apply(scaled_global_count_matrix, 1, max) >= user_inputs()$scale_min, ]
    scaled_global_count_matrix <- scaled_global_count_matrix[!rownames(scaled_global_count_matrix) %in% marker_clash_reactive(), ]
    
    if (nrow(scaled_global_count_matrix) == 0) {
      return(NULL)
    }
    
    dist_matrix <- dist(t(scaled_global_count_matrix), method = "euclidean")
    hc <- hclust(dist_matrix, method = "ward.D2")
    
    predicted_labels <- cutree(hc, k = 6)
    silhouette_values <- silhouette(predicted_labels, dist_matrix)
    avg_silhouette <- mean(silhouette_values[, 3])
    
    matrix_columns = tstrsplit(colnames(subset_count_matrix), "_", fixed=TRUE)[[1]]
    all_markers = rownames(scaled_global_count_matrix)
    
    inverse_gene_cell_hash <- list()
    
    for (key in names(gene_name_lists_by_cell)){
      for (gene in gene_name_lists_by_cell[[key]]){
        inverse_gene_cell_hash[[gene]] <- key
      }
    }
    
    list(
      heatmap_matrix = scaled_global_count_matrix,
      dist_matrix = dist_matrix,
      clustering = hc,
      predicted_labels = predicted_labels,
      avg_silhouette = avg_silhouette,
      subset_count_matrix = subset_count_matrix,
      matrix_columns = matrix_columns,
      all_unique_markers = all_markers,
      inverse_gene_cell_hash = inverse_gene_cell_hash
    )
  })
  
    
    output$heatmapPlot <- renderPlot({
      #data <- heatmap_data()
      
      # if (is.null(data)) {
      #   return(NULL)  # Exit if no valid matrix
      # }
      #marker_lists            <- user_inputs()$gene_name_lists  # Extract marker list reactively
      #print(marker_lists)
      gene_name_lists_by_cell <- user_inputs()$gene_name_lists_by_cell
      # specificity_vector      <- user_inputs()$specificity_vector  # ✅ Access reactively
      # foldchange_vector       <- user_inputs()$foldchange_vector
      # padj_vector             <- user_inputs()$padj_vector
      
      all_markers <- heatmap_data()$all_unique_markers
      
      
      
      
      log2_scaled_global_count_matrix         <- heatmap_data()$heatmap_matrix
      df_row                                  <- data.frame(all_markers)
      rownames(df_row)                        <- all_markers
      colnames(df_row)                        <- "Marker"
      
      inverse_gene_cell_hash <- heatmap_data()$inverse_gene_cell_hash
      
      for(marker in df_row$Marker){
        df_row[marker, "color"] <- paste0(inverse_gene_cell_hash[[marker]], " marker")
      }
      df_row <- df_row[-c(1)]

    
      df_col                                  <- data.frame(heatmap_data()$matrix_columns)
      rownames(df_col)                        <- colnames(heatmap_data()$subset_count_matrix)
      colnames(df_col)                        <- "Cell Type"
      
      # Load RColorBrewer for color palettes
      
      colors <- c("red", "green", "blue", "orange", "violet", "brown")#brewer.pal(6, "Set1")
      
      # Define your colors as a named vector
      color_vector <- c(
        "Astro marker"   = colors[1],
        "Excite marker"  = colors[2],
        "Inhibit marker" = colors[3],
        "Micro marker"   = colors[4],
        "Oligo marker"   = colors[5],
        "OPC marker"     = colors[6]
      )
      
      # Define your colors as a named vector
      cell_vector <- c(
        "Astro"   = colors[1],
        "Excite"  = colors[2],
        "Inhibit" = colors[3],
        "Micro"   = colors[4],
        "Oligo"   = colors[5],
        "OPC"     = colors[6]
      )
      
      
      annotation_colors <- list()
      
      annotation_colors$color       <- color_vector
      annotation_colors$`Cell Type` <- cell_vector
      
      ari <- adjustedRandIndex(heatmap_data()$matrix_columns, heatmap_data()$predicted_labels)
      
      grid.newpage()  # Start a new graphical page
      
      #req(log2_scaled_global_count_matrix, heatmap_data()) 
      
      # print(dim(log2_scaled_global_count_matrix))
      # print(dim(df_col))
      # print(dim(df_row))
      
      pheatmap_obj <- pheatmap(
        log2_scaled_global_count_matrix,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        annotation_col = df_col,
        annotation_row = df_row,
        annotation_colors = annotation_colors,
        fontsize_col = 5,
        fontsize_row = 10,
        fontsize = 25,
        cellwidth = 3,
        cellheight = 9,
        clustering_method = "ward.D2",
        legend = TRUE,
        color = colorRampPalette(c("blue", "white", "red"))(50),
        angle_col = 90
      )
      
      title_grob <- textGrob(
        paste0("Average Silhouette score: ", round(heatmap_data()$avg_silhouette, 2),
               ", ARI score: ", round(ari, 2)), 
        gp = gpar(fontsize = 20, fontface = "bold")
      )
      
      # Arrange title + heatmap in one plot
      grid.arrange(title_grob, pheatmap_obj$gtable, ncol = 1, heights = c(0.1, 1))
      
      
    })
    
    # observe({
    #   print("User Inputs Updated:")
    #   print(studies_thresholded())  # Print to console
    #   #print(str(user_inputs()))
    # })
    
    
    
    # post_hoc_corrected_gene_list_by_cell_type <- list()
    # 
    # for (cell in cell_names){
    #   
    #   post_hoc_corrected_gene_list_by_cell_type[[cell]] <- intersect(gene_name_lists_by_cell[[cell]], 
    #                                                                  rownames(scaled_global_count_matrix))
    #   
    # }
    # 
    # # Display a preview of the data
    # output$dataPreview <- renderTable({
    #   req(dataset())
    #   head(dataset(), 10)
    # })
    
  

    table_data <- reactive({
      
      post_hoc_corrected_gene_list_by_cell_type <- list()
      scaled_global_count_matrix <- heatmap_data()$heatmap_matrix
      gene_name_lists_by_cell    <- user_inputs()$gene_name_lists_by_cell
      
      print(gene_name_lists_by_cell)
      print(rownames(scaled_global_count_matrix))
      
      for (cell in cell_names){

        post_hoc_corrected_gene_list_by_cell_type[[cell]] <- intersect(gene_name_lists_by_cell[[cell]],
                                                                       rownames(scaled_global_count_matrix))

      }
      
      print(post_hoc_corrected_gene_list_by_cell_type)
      
      df <- data.frame(
        Cell_Type = names(post_hoc_corrected_gene_list_by_cell_type),
        Genes = sapply(post_hoc_corrected_gene_list_by_cell_type, paste, collapse = ", "),
        Count = sapply(post_hoc_corrected_gene_list_by_cell_type, length)# Convert to comma-separated strings
      )
      
      return(df)  # Return the dataframe
    })
    
    
    output$gene_table <- renderDT({
      
      datatable(table_data(), 
                options = list(autoWidth = TRUE, scrollX = TRUE), 
                escape = FALSE) %>%
        formatStyle(
          columns = "Genes", # Replace with your column name
          whiteSpace = "normal",
          wordWrap = "break-word"
        )
    })
    
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
