Media independent information handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To handle media independent data access, the :doxy:`GUrl` base class has been 
implemented that represent a media in an abstract way. The class has
abstract ``open``, ``read``, ``write`` and ``close`` method to open
a media, read from it or write to it, and to close the media. 

Support is implemented so far for file and string media, but in the future,
direct access to ressources over the internet may become possible.
A file media is implemented by the :doxy:`GUrlFile` class, while a string
media is implemented by the :doxy:`GUrlString` class.

An example of a class making using of :doxy:`GUrl` is the :doxy:`GXml` class. Look
at the following code:

**Python**

.. code-block:: python
   :linenos:

   xml = gammalib.GXml()
   xml.append(gammalib.GXmlElement('dummy'))
   file  = gammalib.GUrlFile('my_file.xml', 'w')
   chain = gammalib.GUrlString()
   xml.write(file)
   xml.write(chain)

**C++**

.. code-block:: cpp
   :linenos:

   GXml xml;
   xml.append(GXmlElement("dummy"));
   GUrlFile   file("my_file.xml", "w");
   GUrlString chain;
   xml.write(file);
   xml.write(chain);

Line 1 declares a XML object and in line 2 we append a dummy element
to it. In line 3 we now create a file named ``my_file.xml`` for which
we allow write access. In line 4 we allocate a string media. We then
write the XML object first into the file in line 5, and then in the string
in line 6. This illustrates how the :doxy:`GUrl` classes can be used to
redirect the same information to different media. Reading from different
media is analoguous.
