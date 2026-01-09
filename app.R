library(shiny)
library(shinyWidgets)
library(readxl)

validate_excel_file <- function(file_path, 
                                expected_sheets = c("REGIONS", "COMPARISONS", "GROUPS"), 
                                expected_columns = list(REGIONS = c("REGION_GROUP_ID","REGION_GROUP_NAME","REGION_ID","CHROM","START","END"),
                                                        COMPARISONS = c("NUM","GROUP_TREAT"),
                                                        GROUPS = c("GROUP_ID","GROUP_DESC","SAMPLE_ID"))) {
  
  ## sheet REGIONS could be also represent as REGIONS_1, REGIONS_2, ...
  
  # Check if file exists
  if (!file.exists(file_path)) {
    return(list(valid = FALSE, message = "File does not exist"))
  }
  
  # Get sheet names
  sheet_names <- excel_sheets(file_path)
  
  # Check missing sheets
  missing_sheets <- c()
  for (sheet in expected_sheets) {
    if (!any(grepl(sheet, sheet_names))) {
      missing_sheets <- c(missing_sheets, sheet) 
    }
  }
  if (length(missing_sheets) > 0) {
    return(list(valid = FALSE, message = paste("Missing sheets:", paste(missing_sheets, collapse = ", "))))
  }
  
  
  # Check columns in each sheet
  for (sheet in sheet_names) {
    tryCatch({
      df <- read_excel(file_path, sheet = sheet)
      actual_columns <- names(df)
      missing_cols <- if (grepl("REGIONS", sheet)) {
        setdiff(expected_columns[["REGIONS"]], actual_columns)
      } else {
        setdiff(expected_columns[[sheet]], actual_columns)
      }
      
      if (length(missing_cols) > 0) {
        return(list(valid = FALSE, message = paste("Sheet", sheet, "missing columns:", paste(missing_cols, collapse = ", "))))
      }
      
      if (nrow(df) == 0) {
        return(list(valid = FALSE, message = paste0("Empty sheet: ", sheet, " 0 rows")))
      }
    }, error = function(e) {
      return(list(valid = FALSE, message = paste("Error reading sheet", sheet, ":", e$message)))
    })
  }
  
  return(list(valid = TRUE, message = "File validation successful"))
}

# Define UI
ui <- fluidPage(
  
  h2("BoxPlotter for ChIPseq Regions-of-Interest", align = "center"),
  
  # Add some custom CSS for better styling
  tags$head(
    tags$style(HTML("
      .content-wrapper {
        max-width: 800px;
        margin: 0 auto;
        padding: 20px;
      }
      .upload-area {
        border: 2px dashed #ccc;
        border-radius: 10px;
        padding: 40px;
        text-align: center;
        margin: 20px 0;
        background-color: #f9f9f9;
        
        display: flex;
        flex-direction: column;
        align-items: center;
      }
      .upload-area:hover {
        border-color: #999;
        background-color: #f0f0f0;
      }
      .btn-submit {
        background-color: #28a745;
        border-color: #28a745;
        color: white;
        font-weight: bold;
        padding: 10px 30px;
      }
      .btn-submit:hover {
        background-color: #218838;
        border-color: #1e7e34;
      }
      .file-info {
        background-color: #e9ecef;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
      }
      
      .shiny-input-container {
        width: 100%;
        max-width: 400px; /* You can change this value */
      }
    "))
  ),
  
  div(class = "content-wrapper",
      # File upload section
      div(class = "upload-area",
          h4("Upload Configuration File"),
          fileInput("config_file", 
                    label = NULL,
                    accept = c(".xlsx", ".xls"),
                    buttonLabel = "Browse...",
                    placeholder = "Choose XLSX file"),
          p("Select an Excel file (.xlsx or .xls) to use as configuration for the boxplotter job.")
      ),
      
      div(
        style = "margin-top: 15px; font-size: 0.9em; color: #555;",
        # p("The configuration file must be an XLSX file containing at least three sheets:"),
        # tags$ul(
        #   tags$li(strong("GROUPS"), "- defines sample groups"),
        #   tags$li(strong("COMPARISONS"), "- specifies comparisons to make"),
        #   tags$li(strong("REGIONS"), "- contains regions of interest (can be split into multiple sheets like REGIONS_1, REGIONS_2)")
        # ),
        p("Download the example file below to see the required format.")
      ),
      
      # Example file download
      div(
        style = "margin-top: 15px; text-align: center;",
        downloadButton("downloadExample",
                       label = "Download Example Configuration File",
                       icon = icon("download"),
                       class = "btn-primary",
                       style = "width: 100%;")
      ),
      
      # Instructions file download
      div(
        style = "margin-top: 5px; text-align: center;",
        downloadButton("downloadInstructions",
                       label = "Download Instructions",
                       icon = icon("download"),
                       class = "btn-primary",
                       style = "width: 100%;")
      ),
      
      
      # File information display
      conditionalPanel(
        condition = "output.fileUploaded",
        div(class = "file-info",
            h5("Uploaded File Information:"),
            verbatimTextOutput("file_info")
        )
      ),
      
      # Submit button
      div(style = "text-align: center; margin: 30px 0;",
          conditionalPanel(
            condition = "output.fileUploaded",
            actionButton("submit_job", 
                         "Submit Boxplotter Job", 
                         class = "btn btn-lg btn-submit",
                         icon = icon("play"))
          )
      ),
      
      # Job status section
      conditionalPanel(
        condition = "output.showStatus",
        hr(),
        h4("Job Status"),
        verbatimTextOutput("job_status")
      )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive values to store job information
  values <- reactiveValues(
    uploaded_file_path = NULL,
    job_id = NULL,
    job_submitted = FALSE,
    file_validated = FALSE
  )

  # Download example file handler
  output$downloadExample <- downloadHandler(
    filename = function() {
      "example_config.xlsx"  # Name of the downloaded file
    },
    content = function(file) {
      file.copy("www/example.xlsx", file)  # Assumes the file is in the www/ folder
    }
  )
    
  # Download instructions file handler
  output$downloadInstructions <- downloadHandler(
    filename = function() {
      "instructions.pdf"  # Name of the downloaded file
    },
    content = function(file) {
      file.copy("www/instructions.pdf", file)  # Assumes the file is in the www/ folder
    }
  )
  
  # Check if file is uploaded
  output$fileUploaded <- reactive({
    return(values$file_validated)
  })
  
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  
  # Display file information
  output$file_info <- renderText({
    req(input$config_file)
    
    file_info <- input$config_file
    file_size <- round(file_info$size / 1024, 2)  # Convert to KB
    
    paste(
      "Filename:", file_info$name,
      "\nFile size:", file_size, "KB",
      "\nFile type:", tools::file_ext(file_info$name),
      sep = " "
    )
  })
  
  # Handle file upload and copy to permanent location
  observeEvent(input$config_file, {
    req(input$config_file)
    
    ## check the consistence of input file
    if (!is.null(input$config_file)) {
      validation_result <- validate_excel_file(input$config_file$datapath)
      
      if (!validation_result$valid) {
        
        removeNotification("submitting")
        showModal(modalDialog(
          title = tags$h3(icon("exclamation-triangle", style = "color: red;"), " Parsing XLSX config file error"),
          tags$div(
            tags$p("The XLSX file is not valid:"),
            tags$strong("Error message:"), tags$code(validation_result$message)
          ),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        
        # showNotification(
        #   paste("Validation failed:", validation_result$message), 
        #   type = "error"
        # )
        values$file_validated <- FALSE
        return()
      }
    }
    
    values$file_validated <- TRUE
    # Create a permanent directory for uploaded files
    upload_dir <- "uploaded_configs"
    if (!dir.exists(upload_dir)) {
      dir.create(upload_dir, recursive = TRUE)
    }
    
    # Copy file to permanent location with timestamp
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    file_ext <- tools::file_ext(input$config_file$name)
    file_name <- tools::file_path_sans_ext(basename(input$config_file$name))
    permanent_filename <- paste0(file_name,"_", timestamp, ".", file_ext)
    permanent_path <- file.path(upload_dir, permanent_filename)
    
    file.copy(input$config_file$datapath, permanent_path, overwrite = TRUE)
    values$uploaded_file_path <- permanent_path
    #showNotification("File uploaded successfully!", type = "success")
  })
  
  # Submit job when button is clicked
  observeEvent(input$submit_job, {
    req(values$uploaded_file_path)
    
    # Show loading notification
    showNotification("Submitting job...", type = "message", duration = NULL, id = "submitting")
    
    # Construct qsub command
    script_path <- "boxplotter.qsub"  # Adjust path as needed
    config_path <- normalizePath(values$uploaded_file_path)
    qsub_command <- paste("qsub", script_path, config_path)
    
    # Execute qsub command
    tryCatch({
      result <- system(qsub_command, intern = TRUE, ignore.stderr = FALSE)
      
      # Parse job ID from qsub output (typically "Your job XXXXX (jobname) has been submitted")
      job_id_match <- regmatches(result, regexpr("\\d+", result))
      
      if (length(job_id_match) > 0) {
        values$job_id <- job_id_match[1]
        values$job_submitted <- TRUE
        
        # Remove loading notification
        removeNotification("submitting")
        
        # Show success modal with job ID
        showModal(modalDialog(
          title = tags$h3(icon("check-circle", style = "color: green;"), " Job Submitted Successfully!"),
          tags$div(
            tags$p(tags$strong("Job ID:"), tags$code(values$job_id, style = "font-size: 1.2em; color: #007bff;")),
            tags$p(tags$strong("Configuration file:"), basename(values$uploaded_file_path)),
            tags$p(tags$strong("Command executed:"), tags$code(qsub_command)),
            tags$hr(),
            tags$p("You can monitor the job status using:"),
            tags$code(paste("qstat -j", values$job_id)),
            tags$p("Your results will be emailed to you as soon as the job finishes. Feel free to close this tab; your job will continue running in the background.")
          ),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        
      } else {
        # Job submission failed
        removeNotification("submitting")
        showModal(modalDialog(
          title = tags$h3(icon("exclamation-triangle", style = "color: red;"), " Job Submission Failed"),
          tags$div(
            tags$p("The qsub command did not return a valid job ID."),
            tags$p(tags$strong("Command:"), tags$code(qsub_command)),
            tags$p(tags$strong("Output:"), tags$pre(paste(result, collapse = "\n")))
          ),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      }
      
    }, error = function(e) {
      # Handle system command errors
      removeNotification("submitting")
      showModal(modalDialog(
        title = tags$h3(icon("exclamation-triangle", style = "color: red;"), " Error Submitting Job"),
        tags$div(
          tags$p("An error occurred while submitting the job:"),
          tags$pre(as.character(e)),
          tags$p(tags$strong("Command attempted:"), tags$code(qsub_command))
        ),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })
  
  # Show job status section
  output$showStatus <- reactive({
    return(values$job_submitted)
  })
  outputOptions(output, "showStatus", suspendWhenHidden = FALSE)
  
  # Display job status
  output$job_status <- renderText({
    if (values$job_submitted && !is.null(values$job_id)) {
      paste(
        "Job ID:", values$job_id,
        "\nStatus: Submitted",
        "\nConfiguration file:", basename(values$uploaded_file_path),
        "\nSubmission time:", Sys.time(),
        "\n\nTo check job status, run: qstat -j", values$job_id,
        "\nFeel free to close this tab; your job will continue running in the background.",
        sep = " "
      )
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)
