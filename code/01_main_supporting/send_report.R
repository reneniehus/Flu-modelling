send_report <- function(params) {
  
  if (params$send_report) {
    
    library(emayili)
    email <- emayili::envelope()
    Sys.sleep(30) # wait for X seconds, allow time for the html to be rendered
    
    email <- email %>%
      emayili::from(addr = 'rene.niehus@ecdc.europa.eu') %>%
      emayili::to( params$report_recipients ) %>% 
      emayili::subject(subject = 'New flu scenario model run complete') %>%
      #emayili::html('message in html') %>% 
      emayili::attachment(path = './code/03_report/report_overview.html') %>% 
      emayili::attachment(path = './code/03_report/fit_flip.pdf')
    
    smtp <- emayili::server(host = "mailgw.ecdcnet.europa.eu",
                            port = 25,
                            insecure = TRUE,
                            reuse = FALSE)
    smtp(email, verbose = TRUE)
    
  }
  
  
  
  
  return(invisible(NULL))
}

