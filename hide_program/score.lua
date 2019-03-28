#!/usr/bin/env lua
package.cpath = package.cpath .. ";/home/sid/thesis_SidReed/hide_program/?.so"
require 'io'
require 'os'
require 'math'
require 'lpeg'
require 'table'
require 'math'
require 'lfs'
require 'string'
require 'hgt'

-- constants
local chunk_length = 500000

function main(arg)
   -- handle command-line arguments
   if #arg == 0 or arg[1] == '--help' then
      print(usage)
      os.exit(0)
   end
   coff = arg[2] and tonumber(arg[2])
   local species = load_species(arg[1]..'/species.newick')
   local ignore = loadstring('return ' ..
                             ((coff and arg[3]) or
                             ((not coff) and arg[2])
                             or '{}'))()

   for i, edge in ipairs(ignore) do
      uidx, vidx = unpack(edge)
      ignore[i] = {species:node(uidx), species:node(vidx)}
   end

   if coff and coff <= 50 then
      print('ERROR: Bootstrapping cutoff has to be higher than 50%.')
      os.exit(1)
   end

   local genes = {}
   for fname in lfs.dir(arg[1]) do
      if string.match(fname, '[%w_-]*%.newick') then
         if not string.match(fname, 'species%.newick') then
            genes[#genes + 1] = arg[1]..'/'..fname end end end

   local scores = species:blank_scoreboard()
   for i, gene in ipairs(genes) do
      io.stderr:write('processing gene ' .. gene .. ' (' .. i .. ' of ' .. #genes .. ')\n')
      local gene_contrib
      if coff then
         gene_contrib = contrib(species, gene, ignore)
      else
         gene_contrib = single_contrib(species, gene, ignore)
      end
      scores = scores + gene_contrib:expt(2)
   end

   print_folded_scores{scores=scores}
end

usage = [=[
Invoke the script like this: lua score.lua <dir> [cutoff] [ignore]

<dir> should contain a file named species.newick with the species tree
in newick format. All other *.newick files in that directory are taken
to be gene trees. In case bootstrapping is needed, put all versions of
each gene-tree in one file and specify a cutoff percentage. The cutoff
must be higher than 50.

If you want to ignore quartets that are explained by some (directed) edges,
specify [ignore] as follows: "{{u1, v1}, {u2, v2}, ...}" , (do not omit
the quotes!) where ui, vi are the indexes of each edge\'s endpoints.
For undirected edges, specify each edge twice, once in each direction.

The output is printed to stdout, so make sure to redirect it if you want
to save it to a file.

Example:
    lua score.lua foo 80 "{{11, 13}, {13, 11}}" > scores.txt
]=]

local P, R, S, V = lpeg.P, lpeg.R, lpeg.S, lpeg.V
local C, Ct = lpeg.C, lpeg.Ct

local space = (S" \r\n\t")^0

local nwtree = P{
   "tree";
   space = (S" \r\n\t")^0,
   namechar = R("az", "AZ", "09", "._"),
   name = C((V"namechar")^1),
   internal =
      Ct(P"(" * V"node" *
         (P"," * V"space" * V"node")^0 * P")"),
   node = V"internal" + V"name",
   tree = V"node" * P";" * V"space"
}

function nwtree_parse(str)
   return lpeg.match(space * nwtree, str)
end

function choose4(n) return (n*(n-1)*(n-2)*(n-3))/24 end
function qcount(tree) return 3*choose4((tree:node_count()+1)/2) end

function load_species(fname)
   local ifile = io.open(fname)
   local text = ifile:read('*a')
   ifile:close()
   local species = lpeg.match(space * nwtree, text)
   return hgt.species_tree(species) end

function load_genes(fname)
   local ifile = io.open(fname)
   local text = ifile:read('*a')
   ifile:close()
   local genes = lpeg.match(space * lpeg.Ct(nwtree^0), text)
   for i, gene in ipairs(genes) do
      genes[i] = hgt.species_tree(gene) end
   return genes end

function print_folded_scores(arg)
   local hiscores = arg.scores:folded_table()
   local ofile = arg.ofile or io.stdout
   table.sort(hiscores, function (a, b) return a.score > b.score end)
   local count = arg.count or #hiscores
   for i, edge in ipairs(hiscores) do
      if i > count then break end
      if edge.score ~= 0 then
         ofile:write(string.format(
           "%f\t%s --> %s\t(%d%%/%d%%)\n", edge.score,
           tostring(edge.u), tostring(edge.v),
                     math.floor(100*(edge.lscore/edge.score)),
                     math.floor(100*(edge.rscore/edge.score))))
      else
         ofile:write(string.format(
           "%f\t%s --> %s\n", edge.score,
                        tostring(edge.u), tostring(edge.v)))
      end end end

function contrib(species, fname, ignore)
   local raw = species:blank_scoreboard()
   local norm = species:blank_scoreboard()
   local genes = load_genes(fname)
   local treshold = math.ceil(((#genes)*coff)/100)
   for i = 1, qcount(genes[1]), chunk_length do
      raw = raw + species:score(
         species:weighted_quartets(genes, i-1, chunk_length):
         bootstrap(treshold), ignore)
      norm = norm + species:score(genes[1]:
         all_quartets(i-1, chunk_length):translate(genes[1], species), ignore)
   end
   return (raw/norm):keep_max() end

function single_contrib(species, fname, ignore)
   local raw = species:blank_scoreboard()
   local norm = species:blank_scoreboard()
   local gene = load_genes(fname)
   if #gene > 1 then
      print('ERROR: Multiple versions of gene '..
            fname..' , but no cutoff given.')
      os.exit(1)
   end
   gene = gene[1]
   for i = 1, qcount(gene), chunk_length do
      raw = raw + species:score(gene:
         all_quartets(i-1, chunk_length):
         filter_match(gene):
         translate(gene, species), ignore)
      norm = norm + species:score(gene:
         all_quartets(i-1, chunk_length):translate(gene, species), ignore)
   end
   return (raw/norm):keep_max() end

if not inhibit then main(arg) end
