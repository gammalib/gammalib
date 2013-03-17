#! /usr/bin/env python
# ==========================================================================
# Scope
#
#   This script provides an example for creating and handling XML documents
#   with GammaLib
#
# Usage
#   ./xml_howto.py
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
from math import *


# =================== #
# Create XML document #
# =================== #
def create_xml():
    """
    Illustrates the creation of an XML document.
    
    This function creates a simple XML document that illustrates how an
    XML document can be constructed with GammaLib. In this example, the
    document creation is done linearly, in the order in which the nodes
    are written into the document.
    
    The append() method is used to append nodes to the XML document root
    and to child nodes. You will see two flavours of the append() method:
    - one that takes a text string as argument, which serves to define
      the name and the attribtes of the child element
    - one that takes a child node as argument, which serves to append
      either a text leaf using the GXmlText constructor, or other markup
      nodes, such as processing instructions (GXmlPI) or comments
      (GXmlComment)
    """
    # Allocate XML document
    xml = GXml()
    
    # Append first child element to XML document root. The first child
    # element contains a single child which is an empty-element tag.
    level1 = xml.append('element type="Measurement"')
    level2 = level1.append('parameter name="Flux" value="1.0"')

    # Append second child element to XML document root. The second child
    # element contains a single child named list, which contains two
    # childs that each contain a text leaf.
    level1 = xml.append('element')
    level2 = level1.append('list')
    level3 = level2.append('string')
    level4 = level3.append(GXmlText("This is a text"))
    level3 = level2.append('integer')
    level4 = level3.append(GXmlText("17"))

    # Append third child element to XML document root. The third child
    # element contains a single Processing Instruction.
    level1 = xml.append('processing-instruction')
    level2 = level1.append(GXmlPI('<?xml-stylesheet type="text/xsl"?>'))

    # Append forth child element to XML document root. The forth child
    # element contains two comments. The first comment is a single-line
    # comment while the second comment runs over multiple lines.
    level1 = xml.append('comment')
    level2 = level1.append(GXmlComment('This is a comment'))
    level2 = level1.append(GXmlComment('This is a\ntwo line comment'))

    # Return XML document
    return xml


# ======================= #
# Manipulate XML document #
# ======================= #
def manipulate_xml(xml):
    """
    Illustrates the manipulation of an XML document.
    
    This function illustrates the manipulation of an XML document by
    using the remove(), insert(), set() and extend() methods.
    """
    # Remove the third child node containing the processing instruction
    # (Note that the third node has the index 2!).
    xml.remove(2)

    # Insert now at the location where we just removed the processing
    # instruction a new child element named 'replacement'. This new
    # child element will have four levels of child nodes.
    level0 = xml.insert(2, GXmlElement('replacement'))
    level1 = level0.append('level1')
    level2 = level1.append('level2')
    level3 = level2.append('level3')
    level4 = level3.append('level4')

    # Set now a text node at the place of the level4 child node (the
    # level4 child node was sitting in the first slot of the level3
    # node, hence we have to set the text using the index 0)
    level4 = level3.set(0, GXmlText('This text replaces the level 4 child node'))

    # We extend the level2 child node by adding four further child nodes
    # containg a text node, one comment and one processing instruction.
    # In total, the level2 child node now has six childs, four of which
    # are elements, one is a comment and one a processing instruction.
    # 
    # We could do this simply using the append method, but we illustrate
    # here another case: the extension of a child node using child nodes
    # of another element. For this purpose we create first a fresh element
    # 'ext'. To this fresh element we append four child nodes which all
    # will contain a text node. We furthermore append one comment and one
    # processing instruction. Once this is setup we call the extend() method
    # to extend the level2 child node.
    ext      = GXmlElement('level2 extension="yes"')
    element1 = ext.append('level3')
    element2 = ext.append('level3')
    element3 = ext.append('level3')
    element4 = ext.append('level3')
    comment  = ext.append(GXmlComment('This is a comment'))
    pi       = ext.append(GXmlPI('<?xml-stylesheet type="text/xsl"?>'))
    element1.append(GXmlText('This is the first element'))
    element2.append(GXmlText('This is the second element'))
    element3.append(GXmlText('This is the third element'))
    element4.append(GXmlText('This is the forth element'))
    level2.extend(ext)

    # Now we replace the text of the last element in the level2 childs.
    # The emphasis is here on 'element'. The comment and the processing
    # instruction are NOT elements, they are special markup tags. The
    # last 'element' is the forth child node!
    #
    # Instead of explicitly quoting the index of the forth child, we
    # determine the number of elements using the elements() method, and
    # we access the last element by using number-1 as the index in the
    # element() method.
    #
    # The text child is then accessed using the [0] index operator. As
    # there is only a single text node in the last element, we use the
    # index 0.
    number = level2.elements()
    last   = level2.element(number-1)
    text   = last[0]
    text.text('This replaces the text of the last element')

    # Now we change the name of the last element in the level2 childs
    # to 'new3' and we benefit to illustrate also how attributes can
    # be set. The name change is done using the name() method, attributes
    # are manipulated using the attribute() method.
    #
    # There are three calls to attribute(). The first adds a new attribute
    # 'type' with a value of 'none'. As this attribute did not exist before,
    # the attribute() method adds it to the element. The same is done by the
    # second call which add the attribute 'case' with value 'mixed'. The
    # third class modifies the attribute 'type' by giving it a new value of
    # 'info'. This illustrates the logic of the attribute() method: if the
    # attribute name does not exist the attribute will be added, otherwise
    # it will be modified.
    last.name('new3')
    last.attribute('type', 'none')
    last.attribute('case', 'mixed')
    last.attribute('type', 'info')

    # Now we make again a text replacement, but instead of chaning the
    # text of the last element we change the text of the last element with
    # name 'level3'. This is not the last in the list as we just changed
    # the name of the last element to 'new3'. This example illustrates
    # how elements can be accessed by name.
    number = level2.elements('level3')
    text   = level2.element('level3', number-1)[0]
    text.text('This replaces the text of the last "level3" element')

    # Return XML document
    return xml


# ================ #
# XML document I/O #
# ================ #
def xml_io(xml):
    """
    Illustrates input/output for an XML document.
    """
    # Save the XML document into a file
    xml.save("xml_howto.xml")

    # Reload the XML document
    xml.load("xml_howto.xml")
    
    # Return
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
    Simulate photons.
    """
    # Dump header
    print("")
    print("*************************************")
    print("* Example for XML document handling *")
    print("*************************************")
    print("")

    # Create XML document
    xml = create_xml()

    # Show XML document
    print("1. create_xml() output")
    print("======================")
    show_xml(xml)

    # Manipulate XML document
    xml = manipulate_xml(xml)

    # Show XML document
    print("2. manipulate_xml() output")
    print("==========================")
    show_xml(xml)
    
    # Perform I/O
    xml = xml_io(xml)
    
    # Show XML document
    print("3. xml_io() output")
    print("==================")
    show_xml(xml)
    