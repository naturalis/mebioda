Big trees
=========

Recursion
---------

A tree is a hierarchical data structure that is often traversed using recursion. For 
example:

- Assume that a `node` object has a name (`node.name`)
- Assume that a `node` object has zero or more other such nodes as children (`node.children`)

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

How to map a tree to tabular data
----------------------------------

![](phylogeny.png)

For a depth-first traversal, all the children for the focal node need to be reachable, so
a simple table design (table `node` in [tree.db](tree.db)) might be:

| id | name |parent|
|----|------|------|
| 1  | A    | 5    |
| 2  | B    | 5    |
| 3  | C    | 6    |
| 4  | D    | 6    |
| 5  | n1   | 7    |
| 6  | n2   | 7    |
| 7  | n3   | 9    |
| 8  | E    | 9    |
| 9  | n4   |      |

This way, for any focal node, its parent, children and siblings can be queried in 
[SQL](https://en.wikipedia.org/wiki/SQL), e.g. given `n3`, retrieve the parent:

```sql
select p.name from node as n, node as p where n.name = 'n3' and n.parent = p.id
-- n4
```
The children:
```sql
select c.name from node as n, node as c where n.name = 'n3' and n.id = c.parent
-- n1
-- n2
```
Sibling (but not self):
```sql
select s.name from node as n, node as s where n.name = 'n3' and n.parent = s.parent and s.name != 'n3'
-- E
```

Topological queries on tabular trees
------------------------------------

SQL queries quickly become cumbersome, for example when attempting to formulate recursive
traversals. However, there are some shortcuts to implement common queries, for example:

- Get all the descendents of the focal node
- Get all the ancestors
- Get the most recent common ancestor for two nodes

Pre-computing additional columns
--------------------------------

![](phylogeny.png)

These can be implemented using additional, pre-computed columns. Here, `left` is an 
integer that was incremented and assigned to the focal node in pre-order during a 
depth-first traversal, `right` was incremented and assigned post-order:

| id | parent | left | right | name |
|----|--------|------|-------|------|
| 2  | 1      | 1    | 13    | n4   |
| 3  | 2      | 2    | 11    | n3   |
| 4  | 3      | 3    | 6     | n1   |
| 5  | 4      | 4    | 4     | A    |
| 6  | 4      | 5    | 5     | B    |
| 7  | 3      | 7    | 10    | n2   |
| 8  | 7      | 8    | 8     | C    |
| 9  | 7      | 9    | 9     | D    |
| 10 | 2      | 12   | 12    | E    |

Query examples
--------------

Selecting the descendents of `n3`:

```sql
select d.name from node as d, node as n where n.name = 'n3' and d.left > n.left and d.right < n.right;
-- A
-- B
-- n1
-- C
-- D
-- n2
```

Selecting the ancestor(s) of `n3`:

```sql
select a.name from node as a, node as n where n.name = 'n3' and a.left < n.left and a.right > n.right; 
-- n4
```

Selecting the MRCA of `A` and `C`:

```sql
select mrca.name 
	from 
		node as mrca, 
		node as a, 
		node as c 
	where 
		a.name='A' and 
		c.name='C' and 
		m.left < a.left and 
		mrca.right > c.right
	limit 1;
```

Exercise
--------
We are going to figure out which of our crop species are most distant from one another. This means
that, in principe, we have to inspect all pairs - so work together.
- Install a SQLite client if `sqlite3` is unavailable
- Download the database version of the PhytoPhylo tree: https://doi.org/10.6084/m9.figshare.5598631
- Using MRCA queries and the `height` column (distance to root) you should be able to fetch the
  distance between a pair.
