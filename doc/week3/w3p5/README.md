Querying the EoL trait bank
===========================

![](eol.png)

The [Encyclopedia of Life](http://eol.org) is a web project that aims to provide 
global access to knowledge about life on earth. It does this by constructing a web
page for each species and aggregating information from a variety of sources on that
page. As such, the way species are uniquely identified is by their
`page_id`. For example, the `page_id` of the [Sea Otter](https://eol.org/pages/46559130)
is `46559130`, as you can see in the URL of the page.

> Exercise 1: find the `page_id` of your model organism on the EoL website

### Triples and subjects

You will notice that for many species (probably your model organism as well) there is
a tab called _data_ on the species page. This leads to all the trait data that are
available on the EoL website. Each individual trait value has a number of different
bits of information associated with it:

- Where the information originally came from (as of late 2018 this is shown on the right)
- What the name of the trait is
- What the trait value is
- More detailed information when you expand the record by clicking the triangle

For example, a value for the trait `geographic distribution includes` might be `Japan`,
on the basis of a number of observations recorded in the `Japan Species List`. The
basic fact that the geographic distribution of this species includes Japan can be 
represented as a 'triple':

![](triple.png)

A triple is like a machine-readable sentence, and it is called that way because it 
consists of three parts: the subject, the predicate, and the object. For the sentence
'the sky is blue', the subject is 'the sky', the predicate is 'is', and the object
is 'blue'. Likewise, the picture shows a subject (`http://example.org/123`), a predicate
(`ex:age`), and an object (`43`). The triple thus means to say that whatever is identified
by the URL has an age of 43.

> Exercise 2: for any EoL trait record, what is the subject?

### Predicates

On several occasions during the course we've discussed the importance of making it
unambiguously clear what we mean when we use a certain term in our data. For example, 
when you combine different data tables, do columns with the same name automatically 
mean they hold the same data? Are columns with different names necessarily different? 
It is impossible to say unless the terms that we use are anchored on a globally unique
definition. We variously used the word 'ontology' or 'controlled vocabulary' for this.

The EoL trait bank uses this in a big way. Every predicate (so, every term about a
species, such as `geographic distribution includes`) is anchored on an ontology. Many
predicates are basically invented and defined by EoL itself, but others are borrowed from
external ontologies. This is a good idea, because we are only going to be able to
link up on the web of data if there is a shared agreement across different data sources
on what our predicates mean. And to make sure that we uniquely identify the predicate, 
it is given a URL (even if that doesn't necessarily mean that the URL resolves to a web
page! The important part is that it is unique across the entire world, which URLs
automatically are.)

> Exercise 3: for the predicate 'Leaf Complexity', what is its URL? What ontology does it come from?

### Objects

The third part of a triple is the object. In the picture of the triple, the object is
the literal value `43`. It is possible that this is ambiguous, because we don't know
right away what the unit is here. Perhaps it is the case that the definition of the
predicate `ex:age` specifies this, for example by stating that this is age in years
(as opposed to days, millions of years, milliseconds). On other occasions we have 
already seen how this is commonly specified, for example when we collected occurrence
data: the terms [dwc:decimalLatitude](https://terms.tdwg.org/wiki/dwc:decimalLatitude)
and [dwc:decimalLongitude](https://terms.tdwg.org/wiki/dwc:decimalLongitude) are clear
enough. (So, yes, GBIF uses ontologies as well, and so do BoLD, GenBank, and many other
databases with biodiversity data.)

In many cases it becomes unwieldy, or really impossible, to cram information about the
object into the predicate such that the object can remain a literal value. Instead, the
object also becomes a term from an ontology. For example, the territorial extent of 
Japan has varied over the years, so whether a species includes Japan in its geographical
distribution depends on what is meant by that, which EoL clarifies by using
[this globally unique anchor](http://www.geonames.org/1861060).

> Exercise 4: for the trait value 'green' (for example as the object of the triple
> where the predicate is 'leaf color'), what is the definition? What ontology does that
> come from?

<h1>Cypher query</h1>
<form action='https://eol.org/service/cypher'>
  <textarea name='query' id='query' style='clear:all;width:100%' rows='5'>MATCH (n:Trait) RETURN n LIMIT 1;</textarea>
  <input type='submit' />
</form>