import nimibook

var book = newBookFromToc("FrogCap", "book"):
  entry("About", "index")
  section("Tutorials", "tutorials/index"):
    entry("Customize Probeset", "customize")
    entry("Setup Environment", "setup")
    entry("Data Processing", "process")
    entry("Data Alignment", "alignment")
    entry("Troubleshooting", "troubleshoot")
  entry("Download Probes", "download")
  entry("Papers Using FrogCap", "papers")

book.git_repository_url = "https://github.com/chutter/FrogCap-Sequence-Capture"
book.preferred_dark_theme = "coal"
book.favicon_escaped = """<link rel="icon" href="data:image/svg+xml,<svg xmlns=%22http://www.w3.org/2000/svg%22 viewBox=%220 0 100 100%22><text y=%22.9em%22 font-size=%2280%22>&#f52e;</text></svg>">"""

nimibookCli(book)