**_Additional setup steps required for the reproducible research / dynamic documents workshop:_**

1. First, launch RStudio (not R). Click File, New File, R Markdown. This will tell you that you need to install additional packages (knitr, yaml, htmltools, caTools, bitops, and rmarkdown). Click "Yes" to install these.
1. Sign up for a free account at **[RPubs.com](http://rpubs.com/)**.
1. If you want to convert to PDF, you will need to install a **LaTeX** typesetting engine. This differs on Mac and Windows. **Note that this part of the installation may take up to several hours, and isn't strictly required for the workshop.**
    - **Windows LaTeX instructions**:
        1. Download the installer [using this link](http://mirrors.ctan.org/systems/win32/miktex/setup/setup-2.9.5721.exe). It is important to use the full installer, not the basic installer. Run the installer .exe that you downloaded.
        1. Run the installer _twice_, making sure to use the Complete, not Basic, installation:
            1. First, When prompted, select the box to "Download MiKTeX." Select the closest mirror to your location. If you're doing this from Charlottesville, the United States / JMU mirror is likely the closest. This may take a while.
            1. Run the installer again, but this time select "Install" instead of "Download." When prompted _"Install missing packages on-the-fly"_, drag your selection up to "Yes."
    - **Mac LaTeX instructions**:
        1. Download the installer .pkg file [using this link](http://tug.org/cgi-bin/mactex-download/MacTeX.pkg). This is a very large download (>2 gigabytes). It can take a while depending on your network speed.
        1. Run the installer package. 