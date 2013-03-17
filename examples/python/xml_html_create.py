#! /usr/bin/env python
# ==========================================================================
# Scope
#
#   This script creates a simple HTML Web page using GammaLib's xml module
#
# Usage
#   ./xml_html_create.py
#
# -------------------------------------------------------------------------
#
# Copyright (C) 2013 Juergen Knoedlseder
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ==========================================================================
from gammalib import *


# ==================== #
# Create HTML document #
# ==================== #
def create_html():
    """
    Illustrates how GammaLib can be used to create HTML documents.
    
    The HTML document is composed of a main element with name 'html'.
    This elements has two childs: a 'head' child for the header and a
    'body' child for the page information.
    
    The header contains metadata in the 'meta' element and a page
    title defined by the 'title' element.
    
    The body contains arbitrary text. In this example there are two
    lines of text with different font formatting and a GammaLib icon
    that is placed in the top-right corner of the page.
    """
    # Allocate XML document
    xml = GXml()
    
    # Creates HTML base
    html = xml.append('html')

    # Writes header
    header = html.append('head')
    meta   = header.append('meta content="text/html; charset=ISO-8859-1" '
                           'http-equiv="content-type"')
    title  = header.append('title')
    text   = title.append(GXmlText('GammaLib'))
    
    # Writes body
    body   = html.append('body')
    image  = body.append('img style="float: right;" alt="GammaLib" '
                         'src="http://a.fsdn.com/allura/p/gammalib/icon"')
    big    = image.append('big')
    bigger = big.append('big')
    span   = bigger.append('span style="font-family: Arial; font-weight: bold;"')
    text   = span.append(GXmlText('Gammalib'))
    text   = image.append(GXmlText('<br>'))
    span   = image.append('span style="font-family: Arial;"')
    text   = span.append(GXmlText('A versatile toolbox for the scientific analysis '
                          'of astronomical gamma-ray data'))
    
    # Return XML document
    return xml


# ================= #
# Show XML document #
# ================= #
def show_xml(xml):
    """
    Show XML document on the screen.
    """
    # Allocate string URL
    url = GUrlString()
    
    # Write XML document in URL
    xml.write(url)
    
    # Print URL buffer
    print(url.string())
    
    # Return
    return


# ======================== #
# Main routine entry point #
# ======================== #
if __name__ == '__main__':
    """
    This script creates a simple HTML Web page using GammaLib's xml module.
    """
    # Create HTML XML document
    xml = create_html()

    # Show XML document
    show_xml(xml)

    # Save HTML document
    xml.save("gammalib.html")
