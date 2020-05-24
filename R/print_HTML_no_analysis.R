#' Print HTML for users who do not want analysis by Seurat
#' @param parameters Parameters used in copilot_prepros function
#' @param cell_stats Cell statistics
#' @param seq_stats Sequencing statistics
#' @param sample_stats Sample statistics
#' @param dir Directory to the plots
#' @param sample.name Sample name
#' @import R2HTML

print_HTML_no_analysis <- function(parameters, cell_stats, seq_stats, sample_stats, dir, sample.name){
  system(paste0('base64 ', dir, '/UMI_hist_plot.png > ', dir, '/UMI_hist.txt'))
  system(paste0('base64 ', dir, '/Gene_hist_plot.png > ', dir, '/Gene_hist.txt'))
  system(paste0('base64 ', dir, '/bc_rank_plot.png > ', dir, '/barcode_rank.txt'))
  b64_uh <- readChar(paste0(dir, '/UMI_hist.txt'), file.info(paste0(dir, '/UMI_hist.txt'))$size)
  b64_gh <- readChar(paste0(dir, '/Gene_hist.txt'), file.info(paste0(dir, '/Gene_hist.txt'))$size)
  b64_bc <- readChar(paste0(dir, '/barcode_rank.txt'), file.info(paste0(dir, '/barcode_rank.txt'))$size)

  target <- HTMLInitFile(dir, filename=paste0(sample.name, '_summary'))
  HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Calibri">', file=target)
  HTML("<div id='Page1'>", file=target)
  HTML("<div class='topnav'>", file=target)
  HTML(paste0('<a class="notactive" ', 'href="#" ', 'onclick="return show(',"'Page2','Page1'",')"', ";>Analysis</a>"), file=target)
  HTML(paste0('<a class="active" ', 'href="#">Summary</a>'), file=target)
  HTML("</div>", file=target)
  HTML("<div class='title'>", file=target)
  HTML.title(paste0(sample.name,' Summary'), HR=1, file = target)
  HTML.title('Processed by COPILOT', HR=3, file = target)
  HTML("</div>", file = target)

  HTML("<div id='wrapper'>", file=target)

  HTML("<div class='boxed' id='left' align='center'>", file=target)
  HTML('<table style="width:100%">', file=target)
  HTML('<caption>',file=target)
  HTML.title('Parameters', HR=3, file=target)
  HTML('</caption>',file=target)
  HTML(paste('<tr> <td style="color:#444444;">', parameters$para, '</td> <td align="right" style="color:#72BE4B;">', parameters$value, '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  HTML('<table style="width:100%">', file=target)
  HTML('<caption>',file=target)
  HTML.title('Cell Stats', HR=3, file=target)
  HTML('</caption>',file=target)
  HTML(paste('<tr> <td style="color:#444444;">', cell_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', cell_stats$value, '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  HTML('<table style="width:100%">', file=target)
  HTML('<caption>',file=target)
  HTML.title('Sequencing Stats', HR=3, file=target)
  HTML('</caption>',file=target)
  HTML(paste('<tr> <td style="color:#444444;">', seq_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', seq_stats$value, '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  HTML('<table style="width:100%">', file=target)
  HTML('<caption>',file=target)
  HTML.title('Sample Stats', HR=3, file=target)
  HTML('</caption>',file=target)
  HTML(paste('<tr> <td style="color:#444444;">', sample_stats$stat, '</td> <td align="right" style="color:#72BE4B;">', sample_stats$value, '</td> </tr>'), file=target)
  HTML('</table>', file=target)
  HTML("</div>", file = target)

  HTML("<div class='boxed' id='right' align='center'>", file=target)
  HTML("<div class='subtitle'>", file=target)
  HTML.title('UMI Counts Histogram', HR=3, file=target)
  HTML("</div>", file = target)
  HTML(paste0("<img src='data:image/png;base64,", b64_uh, "' width=90%>"), file=target)
  HTML("<div class='subtitle'>", file=target)
  HTML.title('Number of Genes Histogram', HR=3, file=target)
  HTML("</div>", file = target)
  HTML(paste0("<img src='data:image/png;base64,", b64_gh, "' width=90%>"), file=target)
  HTML("<div class='subtitle'>", file=target)
  HTML.title('Barcode Rank Plot', HR=3, file=target)
  HTML("</div>", file = target)
  HTML(paste0("<img src='data:image/png;base64,", b64_bc, "' width=90%>"), file=target)
  HTML("</div>", file = target)

  HTML("</div>", file = target)
  HTML("</div>", file = target)

  HTML('<div id="Page2" style="display:none">', file=target)
  HTML("<div class='topnav'>", file=target)
  HTML(paste0('<a class="active" ', 'href="#">Analysis</a>'), file=target)
  HTML(paste0('<a class="notactive" ', 'href="#" ', 'onclick="return show(',"'Page1','Page2'",')"',";>Summary</a>"), file=target)
  HTML("</div>", file=target)
  HTML("<div class='title'>", file=target)
  HTML.title(paste0(sample.name,' Analysis'), HR=1, file = target)
  HTML.title('Processed by COPILOT', HR=3, file = target)
  HTML("</div>", file = target)

  HTML("<div id='wrapper'>", file=target)
  HTML("<div class='boxed' id='center' align='center'>", file=target)
  HTML("</div>", file = target)
  HTML("</div>", file = target)
  HTML("</div>", file = target)
  HTML('<style type="text/css">
       .title {
       background-color: #003B6D;
       color: white;
       padding: 0px;
       padding-left: 20px;
       padding-bottom: 0px;
       position: fixed;
       top: 0;
       left: 0;
       z-index: 100;
       width: 100%;
       text-align: left;
       LINE-HEIGHT:20px;
       margin: 0px;
       }
       caption {
       text-align: left;
       font-weight: bolder;
       }
       .subtitle{
       #background-color: #ff9900;
       padding-top: 36px;
       margin: 10px;
       text-align: left;
       font-weight: bolder;
       }
       .boxed {
       #background-color: #ff9900;
       #border: 1px solid #868D96;
       padding: 10px;
       margin: 30px;
       }
       h1 {
       font-family: Calibri;
       font-size: 40px;
       }
       h2 {
       font-family: Calibri;
       font-size: 32px;
       }
       h3 {
       font-family: Calibri;
       font-size: 26px;
       }
       #wrapper {
       display: flex;
       align-items: flex-start;
}
#left {
width: 50%;
}
#right {
width: 50%;
}
#center {
width: 100%;
}
table {
#background-color: #ff4040;
padding: 0px;
margin: 0px;
font-family: "Calibri";
font-size: 22px;
#border: 1px solid #CDCDCD;
}
td {
border-bottom: 1px solid #CDCDCD;
}
#mathplayer{
height: 80px;
}
/* Add a black background color to the top navigation */
.topnav {
background-color: transparent;
overflow: hidden;
top: 56px;
right: 0;
width: 100%;
position: fixed;
padding: 0px;
margin: 0px;
z-index: 999;
text-decoration: none;
}

/* Style the links inside the navigation bar */
.topnav a {
#border-bottom: 1px solid white;
float: right;
color: white;
text-align: center;
text-decoration: none;
font-size: 26px;
font-family: Calibri;
padding: 20px;
padding-top: 10px;
padding-bottom: 6px;
margin: 0px;
}

/* Change the color of links on hover */
.topnav a:hover {
#border-bottom: 1px solid black;
background-color: white;
color: black;
float: right;
color: "#EBEDF3";
text-align: center;
text-decoration: none;
font-size: 26px;
font-family: Calibri;
padding: 20px;
padding-top: 10px;
padding-bottom: 6px;
margin: 0px;
}

/* Add a color to the active/current link */
.topnav a.active {
#border-bottom: 1px solid white;
background-color: white;
color: black;
float: right;
color: "#EBEDF3";
text-align: center;
text-decoration: none;
font-size: 26px;
font-family: Calibri;
padding: 20px;
padding-top: 10px;
padding-bottom: 6px;
margin: 0px;
}
</style> </head>', file=target)
  HTMLEndFile()
  system(paste0("sed -i '/<hr/d' ", dir, "/", sample.name, "_summary.html"))
  system(paste0("sed -i '/<font/d' ", dir, "/", sample.name, "_summary.html"))
  system(paste0("sed -i '/R2HTML/d' ", dir, "/", sample.name, "_summary.html"))
  system(paste0("sed -i 's/nequations=0;/function show(shown, hidden) {document.getElementById(shown).style.display=",'"block"',';document.getElementById(hidden).style.display="none"',";return false;}/g' ", dir, "/", sample.name, "_summary.html"))
  }
