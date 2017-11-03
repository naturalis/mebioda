Big trees
=========

Recursion
---------

A tree is a hierarchical data structure that is often traversed using recursion. For 
example:

```python
# Start at the root
depth_first(root)

# Function definition. This calls
# itself for every child of the input node
def depth_first(node):
	
	# pre-order
	print 'PRE: ', node.name
	
	# recurse further
	for child in node.children:
		depth_first(child)
	
	# post-order
	print 'POST: ', node.name

```

Depth-first traversal
---------------------

![](phylogeny.png)

A [depth-first](https://en.wikipedia.org/wiki/Depth-first_search) traversal explores as
deeply as possible before back tracking, and so the  result of the pre- and post-order 
actions would be (with this [script](recursion.pl) and this [input file](tree.nex)):

```
$ perl recursion.pl 
PRE: n4
PRE: n3
PRE: n1
PRE: A
PRE: B
POST: B
POST: A
PRE: n2
PRE: C
PRE: D
POST: D
POST: C
POST: n2
POST: n1
PRE: E
POST: E
POST: n3
POST: n4
```

Breadth-first traversal
-----------------------

![](phylogeny.png)

A breadth-first algorithm:

```python
def breadth_first(node):
	print 'PRE: ', node.name
	
	# first explore siblings recursively
	if node.sibling:
		breadth_first(node.sibling)
	
	# then move deeper
	if node.first_child:
		breadth_first(node.first_child)
	
	print 'POST: ', node.name
```

Would produce the following:

```
PRE: n4
PRE: n3
PRE: E
POST: E
PRE: n1
PRE: n2
PRE: C
PRE: D
POST: D
POST: C
POST: n2
PRE: A
PRE: B
POST: B
POST: A
POST: n1
POST: n3
POST: n4
```