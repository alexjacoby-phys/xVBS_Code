#########################################################################
#
#  Density Matrix Renormalization Group (and other methods) in julia (DMRjulia)
#                              v0.1
#
#########################################################################
# Made by Thomas E. Baker (2018)
# See accompanying license with this program
# This code is native to the julia programming language (v1.1.0) or (v1.5)
#

"""
    Module: QN

Base definitions for quantum numbers
"""
module QN
using ..tensor

  """
      Qnum

  Abstract types for quantum numbers
  """
  abstract type Qnum end

  """
      fermionQnum

  Abstract types for fermionic quantum numbers that includes parity
  """
  abstract type fermionQnum <: Qnum end
  export Qnum,fermionQnum


  import Base.>
  import Base.<
  import Base.>=
  import Base.==
  import Base.!=
  import Base.<=
  import Base.+
  import Base.inv
  import Base.isequal
  import Base.isless
  import Base.isgreater
  ##All comparison operator are available for every derived Qnum Types
  ## if the operator > and the operator == are available.
  @inline function !=(a::T,b::T)::Bool where {T<:Qnum}
    return !(a==b)
  end
  @inline function <(a::T,b::T)::Bool where {T<:Qnum}
    return a!=b && !(a>b)
  end
  @inline function >=(a::T,b::T)::Bool where {T<:Qnum}
    return !(a<b)
  end
  @inline function <=(a::T,b::T)::Bool where {T<:Qnum}
    return !(a>b)
  end

  @inline function isless(a::T,b::T)::Bool where {T<:Qnum}
    return a<b
  end
  @inline function isgreater(a::T,b::T)::Bool where {T<:Qnum}
    return a>b
  end
  @inline function isequal(a::T,b::T)::Bool where {T<:Qnum}
    return a==b
  end

  @inline function parity(a::T) where {T <: fermionQnum}
	return a.parity.val&1
  end
  export parity

  """
      U1(q)

  Generate a `U1` quantum number with value `q` (default 0)

  #Example:

  ```julia
  U1() #zero quantum number

  U1(2)
  ```

  #Note:

  Encode U1 quantum numbers: an infinite discrete abelian group
  The constructor with no argument creates the neutral element.
  The group operation is the addition (+)
  For any given element the inverse can be computed with inv
  or the element can be inverted with inv!
  By definition "inv(a) + a" is the neutral element for any element "a" of the group.
  Additionnally comparison operator <,>,>=,<=,==,!= are defined. The comparators are not meaningful
  from a group theoretical point of view.

  See also: [`U1`](@ref) [`@composeQNs`](@ref)
  """
  mutable struct U1<:Qnum
    val::Int64
  end

  @inline function U1()::U1
    return U1(0)
  end

  @inline function U1(a::U1)::U1
    return U1(a.val)
  end
  export U1

  @inline function >(a::U1,b::U1)::Bool
    return a.val>b.val
  end
  @inline function <(a::U1,b::U1)::Bool
    return a.val<b.val
  end
  @inline function ==(a::U1,b::U1)::Bool
    return a.val==b.val
  end

  import Base.copy
  @inline function copy(a::U1)::U1
    return U1(a.val)
  end

  @inline Base.@pure function +(a::U1,b::U1)::U1
    ap = copy(a)
    return add!(ap,b)
#    return U1(a.val+b.val)
  end

  import .tensor.add!
  @inline Base.@pure function add!(a::U1,b::U1)::U1
    a.val += b.val
    a
    # nothing
  end



  @inline function inv(a::U1)::U1
    return U1(-a.val)
  end
  export inv

  @inline function inv!(a::U1)::U1
	a.val *=-1;
	a
#   nothing
  end


  @inline function hash(a::U1)::number
    return a.val
  end
  export inv!

  """
      Zn{N}(q)

  Generate a `Zn` quantum number of modulus `N` with value `q` (default 0)

  #Example:

  ```julia
  Zn{2}() #zero quantum number

  Zn{4}(2)
  ```

  # Use:
  Encode Z_n quantum numbers: a finite abelian group with n elements
  The constructor with no argument creates the neutral element.
  The group operation is the addition (+)
  For any given element the inverse can be computed with inv
  or the element can be inverted with inv!
  By definition "inv(a) + a" is the neutral element for any element "a" of the group.
  additionnally comparison operator <,>,>=,<=,==,!= are defined. The comparators are not meaningful
  from a group theoretical point of view.

  See also: [`U1`](@ref) [`@composeQNs`](@ref)
  """
  mutable struct Zn{T} <: Qnum
    val::UInt64

    function Zn{T}(a) where {T}
#      @assert(typeof(T) <: Integer, "The template parameter must be an integral value")
      return a < T ? new{T}(a) : error("Value in Zn{$T} must be smaller than $T")
    end
  end

  function Zn{N}(a::Zn{N})::Zn{N} where{N}
    return Zn{N}(a.val)
  end
  function Zn{N}()::Zn{N} where {N}
    return Zn{N}(0)
  end
  export Zn

  @inline function >(a::Zn{N},b::Zn{N})::Bool where{N}
    return a.val>b.val
  end
  @inline function <(a::Zn{N},b::Zn{N})::Bool where{N}
    return a.val<b.val
  end
  @inline function ==(a::Zn{N},b::Zn{N})::Bool where{N}
    return a.val==b.val
  end

  @inline Base.@pure function add!(a::Zn{N},b::Zn{N})::Zn{N} where{N}
    a.val += b.val
	# a.val %= N
	# this conditionnal is less costly than the full modulo.
	# It is correct if both b.val and a.val are smaller than N
    if a.val >= N
      a.val -= N
    end
    a
    # nothing
  end

  @inline function copy(a::Zn{N})::Zn{N} where{N}
    return Zn{N}(a.val)
  end

  @inline Base.@pure function +(a::Zn{N},b::Zn{N})::Zn{N} where{N}
    #implemented in term of the in place addition.
    ap = copy(a)
    return add!(ap,b)
  end
  export add!

  @inline function inv(a::Zn{N})::Zn{N} where{N}
    if a.val == 0 #The neutral element is always its own inverse.
      return Zn{N}()
    end
    return Zn{N}(N-a.val)
  end
  @inline function inv!(a::Zn{N})::Zn{N}  where{N}
    if a.val != 0 #The neutral element is always its own inverse.
      a.val = (N-a.val);
	end
	a
#   nothing
  end

  @inline function hash(a::Zn{N}) where{N}
    return a.val
  end

  #no varargs template in julia, for some reason... so the simple solution is not possible.
  #must resort to macro...
  # struct compQn{T...}
  # 	struct wrap{W}
  # 		val::W
  # 	end
  # 	TT = (wrap{t} for t in T)
  # 	Content::tuple{TT...}
  # end

  """
      @composeQNs("name",types...)

  creates a quantum number with `name` (String) with `types` in order; pre-defined types: U1, Zn, parity

  # Example:

     @composeQNs("spin",U1) #spins

     @composeQNs("fermion",U1,U1) #fermions

  # Use:

    macro to automatically compose quantum numbers.
      The macro composeQNs_printcode prints the generated code for inspection if you want
      to look at the code this macro would add to yours.
    Supply to the macro a name for your new quantum number type and
    all the quantum number that must be composed to make it.
    The composite quantum number has all the same functions available as a simple or handwritten quantum,
    like all type derived from Qnum should: the group operation: +, copy inverse: inv,in-place inverse: inv!,
    and the comparison operators: ==, !=, >, <, >=, <=.
	suppose we have a system where the number of particle, total spin, and some abstract Z2 quantity are conserved.
	the number of particle N is a U1 number, the total spin S is a U1 number.
    The composite quantum number can be created by calling:
      @composeQNs("U1xU1xZ2",U1,U1,Zn{2})
    after having called this macro, a new type named "U1xU1xZ2" is available.
    we can create such a quantum number with total particle 4, spin 0, and Z2 0 with the call
      N4S0P0 = U1xU1xZ2(4,0,0)
    Be sure to use a name that is unique to your current scope to avoid name clashes.
    Name all your composed QNs with a unique name.

  See also: [`@composeQNs`](@ref) [`U1`](@ref) [`Zn`](@ref) [`parity`](@ref) [`@composeQNs_printcode`](@ref) [`@fermion_composeQNs`](@ref) [`_generatecompositestring`](@ref)
  """
  macro composeQNs(name::String,Types...)
    # _generatecompositestring
    output = _generatecompositestring(name,Types...)
    EExp = Meta.parse(output)
    return eval(EExp)
  end
  export @composeQNs

  """
      @fermion_composeQNs("name",types...)

  creates a quantum number with `name` (String) with `types` in order; pre-defined types: U1, Zn, parity

  # Use:

	 macro to create fermionic quantum numbers.
	 this macro also has an extra quantum number for the parity of the wave function. It is the first value stored within.
	 if we have an object A<:fermionQnum, its parity value can obtained by calling parity(A).
	 Just like with the regular composeQNs macro, you create a new composite quantum number type by calling
		@fermion_composeQNs(<name>,<QNs to compose>...)
	 to create a composite quantum number type with the name of your choice.
   Using this sort of quantum number make it possible to apply fermionic field operator on an qMPS state completly locally.
   
  See also: [`@composeQNs`](@ref) [`U1`](@ref) [`Zn`](@ref) [`parity`](@ref) [`@composeQNs_printcode`](@ref) [`_generatecompositestring`](@ref)
  """
  macro fermion_composeQNs(name::String,Types...)
	output = _generatecompositestring(name,Types...,Inheritfrom=fermionQnum,valfield=["parity"])
	EExp = Meta.parse(output)
    return eval(EExp)
  end
  export @fermion_composeQNs

  """
      @composeQNs_printcode("name",types...)

  prints code created when defining a quantum number type

  # Use:

  macro to create fermionic quantum numbers.
  this macro also has an extra quantum number for the parity of the wave function. It is the first value stored within.
  if we have an object A<:fermionQnum, its parity value can obtained by calling parity(A).
  Just like with the regular composeQNs macro, you create a new composite quantum number type by calling
  @fermion_composeQNs(<name>,<QNs to compose>...)
  to create a composite quantum number type with the name of your choice.
  Using this sort of quantum number make it possible to apply fermionic field operator on an qMPS state completly locally.

  See also: [`@composeQNs`](@ref) [`U1`](@ref) [`Zn`](@ref) [`parity`](@ref) [`@fermion_composeQNs`](@ref) [`_generatecompositestring`](@ref)
  """
  macro composeQNs_printcode(name::String,Types...)
    output = _generatecompositestring(name,Types...)
    println(output)
    nothing
  end
  export @composeQNs_printcode


  """
      _generatecompositestring(name,Types[,Inheritfrom=,valfield=])

  Helper function to make quantum numbers with `@composeQNs`, `@composeQNs` and `composeQNs_printcode`; defines basic operations =,>,+, etc.

  See also: [`@composeQNs`](@ref) [`U1`](@ref) [`Zn`](@ref) [`parity`](@ref) [`@@composeQNs_printcode`](@ref) [`@fermion_composeQNs`](@ref)
  """
  @inline function _generatecompositestring(name,Types...;Inheritfrom = Qnum,valfield = String[])
	if length(valfield) < length(Types)
		lvf = length(valfield)
		for (i,T) in enumerate(Types[lvf+1:end])
			push!(valfield,"val$(i+lvf)")
		end
	elseif length(valfield) > length(Types)
		error("more names than there are quantum number to compose")
	end
#	@assert(Inheritfrom <: Qnum)
    structstr = "mutable struct $name <:$Inheritfrom\n"
    funadd = "@inline function +(a::$name,b::$name)::$name\n	return $name("
    fungreater = "@inline function >(a::$name,b::$name)::Bool\n\t"
    funequal = "@inline function ==(a::$name,b::$name)::Bool\n"
    funinv = "@inline function inv(a::$name)::$name\n	return $name("
    funinv! ="@inline function inv!(a::$name)::$name\n"
    funadd! ="@inline function add!(a::$name,b::$name)::$name\n"
    funNeut = "@inline function $name()\n	return $name("
    funpropercopy = "@inline function copy(a::$name)\n        $name("
    funcopy = "@inline function copy!(a::$name,b::$name)\n for w = 1:$(length(valfield))\n"
#    funcinv = "@inline function copyinv!(a::$name,b::$name)::$name\n"
#    funaddinv = "@inline function addinv!(a::$name,b::$name)::$name\n"
    funsize = "@inline function length(Q::$name)\n\t return $(length(valfield))\nend"
    fungetindex = "@inline function getindex(a::$name,i::Integer)\n        if i == 1"
     funsetindex = "@inline function setindex!(a::$name,b::Integer,i::Integer)\n        if i == 1"
    constructname = "$name("
    constructbody = "new("
    copyconst = "@inline function $name(a::$name)\n\t$name("
    counter = 2

    funcopy = "$(funcopy) \ta[w] = b[w]\n"
    funequal = "$funequal	for w = 1:$(length(valfield)) \n if a[w] != b[w]\n\t\treturn false\n\tend\n end\n"
    for (field,T) in zip(valfield,Types)
      # @assert $T <: Qnum
      structstr = "$structstr	$field::$T\n"
      funadd = "$(funadd)a.$field+b.$field, "
      fungreater = "$(fungreater)if a.$field != b.$field \n\t\treturn a.$field>b.$field\n\tend\n\t"
      funinv = "$(funinv)inv(a.$field), "
      funinv! = "$funinv!	inv!(a.$field);\n"
      funadd! = "$funadd!	add!(a.$field,b.$field);\n"
      funNeut = "$(funNeut)$T(), "
      funpropercopy = "$(funpropercopy) a.$field,"
#      funcinv = "$(funcinv) \ta.$field = inv(b.$field)\n"
#      funaddinv = "$(funaddinv) \tadd!(a.$field, inv(b.$field))\n"
  #		funsize = "@inline function size(Q::$name,i::Integer...)\n"
      fungetindex = "$(fungetindex)\n\treturn a.$field.val\n\telseif i == $counter "
       funsetindex = "$(funsetindex)\n\t a.$field.val = b\n\telseif i == $counter "
      constructname = "$(constructname)l$field, "
      constructbody = "$(constructbody)$T(l$field), "
      copyconst = "$(copyconst) a.$field, "
      counter += 1
    end
    constructor = "\t$(constructname[1:end-2])) = $(constructbody[1:end-2]))"
    structstr = "$(structstr)\n$(constructor)\nend\n"
    funadd = "$(funadd[1:end-2]))\nend"
    fungreater = "$(fungreater)\n	return false\n end"
    funequal = "$funequal	return true\n end"
    funinv = "$(funinv[1:end-2]))\nend"
    funinv! = "$funinv!\n\ta\nend"
    funadd! = "$funadd!\n\ta\nend"
    funNeut = "$(funNeut[1:end-2]))\nend"
    funpropercopy = "$(funpropercopy[1:end-1]) )\nend\n"
    funcopy = "$(funcopy) \n end\n\ta\nend\n"
    fungetindex = "$(fungetindex[1:end-14])end\nend"
     funsetindex = "$(funsetindex[1:end-14])end\nend"
    copyconst = "$(copyconst[1:end-2]))\nend"
#    funcinv = "$(funcinv) \n\ta\nend\n"
#    funaddinv = "$(funaddinv) \n\ta\nend\n"
    # return structstr,funadd,fungreater,funequal,funinv,funinv!,funNeut
    output = """begin
    $(structstr)\n
    $(copyconst)\n
    import Base.>;import Base.<;import Base.==;import Base.+;import Base.length;import Base.getindex;import Base.setindex!; import Base.copy; import Base.copy!\n
    $funadd\n
    $fungreater\n
    $funequal\n
    $funinv\n
    $funinv!\n
    $funadd!\n
    $funNeut\n
    $funsize\n
    $fungetindex\n
    $funsetindex\n
    $funcopy\n
    $copyconst\n
    $funpropercopy\n
    export $name\n
    end"""
    return output
  end

end
