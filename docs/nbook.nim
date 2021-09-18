import nimibook

var book = newBookFromToc("FrogCap", "book"):
  entry("Introduction", "index")
  section("Tutorials", "tutorials/index"):
    entry("Customize Probeset", "customize")
    entry("Setup Environment", "setup")
    entry("Data Processing", "process")
    entry("Data Alignment", "alignment")
    entry("Troubleshooting", "troubleshoot")
  entry("Download Probes", "download")
  entry("Papers Using FrogCap", "papers")

book.git_repository_url = "https://github.com/chutter/FrogCap-Sequence-Capture"

nimibookCli(book)