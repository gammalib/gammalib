.. _sec_xml:

XML file interface
------------------

Overview
~~~~~~~~

The followinf figure presents an overview over the C++ classes of the XML
module and their relations.

.. _fig_uml_xml:

.. figure:: uml_xml.png
   :width: 80%
   :align: center

   *XML module*

The XML module provides classes that allow creation, writing and reading of
files in the Extensible Markup Language (XML) format.
The central class of the XML module is the abstract :doxy:`GXmlNode` base
class which represent a node of the XML document. The essential property
of :doxy:`GXmlNode` is that it may contain a list of :doxy:`GXmlNode` objects,
allowing creating a complex tree of data structures. Nodes in a XML
document may be either elements, comments, some text of a processing
instruction. These different node types are implemented by the
:doxy:`GXmlElement`, :doxy:`GXmlComment`, :doxy:`GXmlText` and :doxy:`GXmlPI` classes,
respectively. A special node if the :doxy:`GXmlDocument` node which provides
the root node of a XML document. This node must exist only once in a XML
tree. The XML file is implement by the :doxy:`GXml` class which contains one
instance of :doxy:`GXmlDocument`.


Creating a XML file
~~~~~~~~~~~~~~~~~~~

The following example illustrates the creation of a XML file
(see ``examples/cpp/createxml/createxml.cpp`` for the source code; the line numbers are
for reference and are not part of the source code):

.. code-block:: cpp
   :linenos:

   GXml xml;
   GXmlComment comment("Now 2 nodes with spatial and spectral info");
   GXmlElement spatial("spatial type=\"Position\"");
   GXmlElement spectral("spectrum type=\"PowerLaw\"");
   GXmlElement text("text");
   GXmlPI      pi("<?process now?>");
   spatial.append(GXmlElement("parameter ra=\"83.0\""));
   spatial.append(GXmlElement("parameter dec=\"22.0\""));
   spectral.append(GXmlElement("parameter prefactor=\"1e-7\""));
   spectral.append(GXmlElement("parameter index=\"-2.1\""));
   text.append(GXmlText("Finish with text"));
   xml.append(comment);
   xml.append(spatial);
   xml.append(spectral);
   xml.append(text);
   xml.append(pi);
   xml.save("my_xml_file.xml");

Below the content of the XML file that will be created by this code:

.. code-block:: xml

   <?xml version="1.0" encoding="UTF-8" standalone="no"?>
   <!--Now 2 nodes with spatial and spectral info-->
   <spatial type="Position">
    <parameter ra="83.0" />
    <parameter dec="22.0" />
   </spatial>
   <spectrum type="PowerLaw">
    <parameter prefactor="1e-7" />
    <parameter index="-2.1" />
   </spectrum>
   <text>Finish with text</text>
   <?process now?>

In line 1 a XML object if allocated. In lines 2-6, one comment node, three
element nodes and one processing instruction node are created. For the
comment node, the comment text is provided in the :doxy:`GXmlComment`
constructor. For the element nodes, the element tag as well as any
attributes (separated by whitespace characters) are provided in the
:doxy:`GXmlElement` constructor. For the processing instruction node,
the instruction including the brackets are provided in the :doxy:`GXmlPI`
constructor. In lines 7-8, two element nodes providing spatial
parameters are appended to the ``spatial`` node. In lines 9-10, parameters
are appended to the ``spectral`` node. And in line 11, a text node is
appended to the ``text`` node.
Finally, the nodes are appended in lines 12-16 to the document root, and
the XML file is saved in line 17.

Note that in the above example the entire XML tree has been constructed
before the nodes were appended to the document root. The reason behind
this approach is that the :doxy:`GXml::append` method creates deep copies of the
nodes provided in the argument, hence manipulation of the node once 
appended requires to retrieve the pointers to the deep copies in the XML
document. The following example illustrates how this can be done:

.. code-block:: cpp

    xml.append(GXmlElement("source type=\"PointSource\""));
    xml.element("spatial", 0)->append(GXmlElement("parameter ra=\"83.0\""));
    xml.element("spatial", 0)->append(GXmlElement("parameter dec=\"22.0\""));

The ``xml.element("spatial", 0)`` method returns a pointer to the first
:doxy:`GXmlElement` node with tag ``spatial`` in the XML document. Now that
we have a pointer to the nodes, elements can be appended to the XML
document using the :doxy:`GXml::append` method.

Alternatively, one can also retrieve the node pointer when the node is 
appended to the XML document:

.. code-block:: cpp

    GXmlNode* node = xml.append(GXmlElement("source type=\"PointSource\""));
    node->append(GXmlElement("parameter ra=\"83.0\""));
    node->append(GXmlElement("parameter dec=\"22.0\""));
 
The :doxy:`GXml::append` method returns in fact the pointer to the deep copy of the
element that has been appended. This pointer can then be used to manipulate
directly the nodes in the XML document.
