#!/usr/bin/env perl

#$Header: /usr/local/ollincvs/Codes/OllinSphere-BiB/prl/arrays.pl,v 1.53 2021/06/18 19:06:05 malcubi Exp $

# This perl script creates the subroutines:
#
# arrays.f90
# accumulate.f90
# allocatearrays.f90
# grabarray.f90
# saveold.f90
# simpleboundary.f90
# symmetries.f90
# syncall.f90
# update.f90
#
# Plus the include files:
#
# boundinterp.inc
# restrict_copy.inc
# restrict_send.inc
# restrict_recv.inc

print "CREATING FILES FOR DEALING WITH ARRAYS ...\n\n";

# Open outfiles.

open(FILE_ARRAYS,">src/auto/arrays.f90") or die "Can't open arrays.f90: $!";
open(FILE_ACCUMULATE,">src/auto/accumulate.f90") or die "Can't open accumulate.f90: $!";
open(FILE_ALLOCATEARRAYS,">src/auto/allocatearrays.f90") or die "Can't open allocatearrays.f90: $!";
open(FILE_GRABARRAY,">src/auto/grabarray.f90") or die "Can't open grabarrays.f90: $!";
open(FILE_SAVEOLD,">src/auto/saveold.f90") or die "Can't open saveold.f90: $!";
open(FILE_SIMPLEBOUNDARY,">src/auto/simpleboundary.f90") or die "Can't open simpleboundary.f90: $!";
open(FILE_SYMMETRIES,">src/auto/symmetries.f90") or die "Can't open symmetries.f90: $!";
open(FILE_SYNCALL,">src/auto/syncall.f90") or die "Can't open syncall.f90: $!";
open(FILE_UPDATE,">src/auto/update.f90") or die "Can't open update.f90: $!";

open(FILE_BOUNDINTERP,">src/auto/boundinterp.inc") or die "Can't open boundinterp.inc: $!";
open(FILE_RESTRICTCOPY,">src/auto/restrict_copy.inc") or die "Can't open restrict_copy.inc: $!";
open(FILE_RESTRICTSEND,">src/auto/restrict_send.inc") or die "Can't open restrict_send.inc: $!";
open(FILE_RESTRICTRECV,">src/auto/restrict_recv.inc") or die "Can't open restrict_recv.inc: $!";
open(FILE_RESTRICTCOSMO,">src/auto/restrict_cosmo.inc") or die "Can't open restrict_cosmo.inc: $!";

# Write beginning of file arrays.f90

print FILE_ARRAYS "! Automatically generated file.  Do not edit!\n\n";
print FILE_ARRAYS "  module arrays\n\n";
print FILE_ARRAYS "  integer, allocatable, dimension(:) :: s\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: t\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: t1\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: t2\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: dr\n";
print FILE_ARRAYS "  real(8), allocatable, dimension(:) :: dt\n\n";
print FILE_ARRAYS "  real(8), dimension (:,:,:), pointer :: grabvar_bound\n\n";
print FILE_ARRAYS "  character(50), allocatable, dimension(:) :: outvars0Darray\n";
print FILE_ARRAYS "  character(50), allocatable, dimension(:) :: outvars1Darray\n";

# Write beginning of file accumulate.f90

print FILE_ACCUMULATE "! Automatically generated file.  Do not edit!\n\n";
print FILE_ACCUMULATE "  subroutine accumulate(l,k,niter,w)\n\n";
print FILE_ACCUMULATE "  use param\n";
print FILE_ACCUMULATE "  use arrays\n\n";
print FILE_ACCUMULATE "  implicit none\n\n";
print FILE_ACCUMULATE "  logical contains\n\n";
print FILE_ACCUMULATE "  integer l,k,niter\n\n";
print FILE_ACCUMULATE "  real(8) w\n\n";

# Write beginning of file allocatearrays.f90

print FILE_ALLOCATEARRAYS "! Automatically generated file.  Do not edit!\n\n";
print FILE_ALLOCATEARRAYS "  subroutine allocatearrays(status)\n\n";
print FILE_ALLOCATEARRAYS "  use param\n";
print FILE_ALLOCATEARRAYS "  use arrays\n";
print FILE_ALLOCATEARRAYS "  use procinfo\n\n";
print FILE_ALLOCATEARRAYS "  implicit none\n";
print FILE_ALLOCATEARRAYS "  logical contains\n";
print FILE_ALLOCATEARRAYS "  character(len=*) status\n\n";
print FILE_ALLOCATEARRAYS "  if (trim(status)=='on') then\n";
print FILE_ALLOCATEARRAYS "     allocate(s(0:Nl-1),t(0:Nl-1),t1(0:Nl-1),t2(0:Nl-1),dr(0:Nl-1),dt(0:Nl-1))\n";
print FILE_ALLOCATEARRAYS "     s = 0; t = 0.d0; dr=0.d0; dt=0d0\n";
print FILE_ALLOCATEARRAYS "  else\n";
print FILE_ALLOCATEARRAYS "     deallocate(t,dr,dt)\n";
print FILE_ALLOCATEARRAYS "  end if\n\n";

# Write beginning of file grabarray.f90

print FILE_GRABARRAY "! Automatically generated file.  Do not edit!\n\n";
print FILE_GRABARRAY "  subroutine grabarray(varname)\n\n";
print FILE_GRABARRAY "  use param\n";
print FILE_GRABARRAY "  use arrays\n";
print FILE_GRABARRAY "  use procinfo\n\n";
print FILE_GRABARRAY "  implicit none\n";
print FILE_GRABARRAY "  logical exists\n";
print FILE_GRABARRAY "  character(len=*) varname\n\n";
print FILE_GRABARRAY "  exists = .false.\n\n";

# Write beginning of file saveold.f90

print FILE_SAVEOLD "! Automatically generated file.  Do not edit!\n\n";
print FILE_SAVEOLD "  subroutine saveold(l)\n\n";
print FILE_SAVEOLD "  use param\n";
print FILE_SAVEOLD "  use arrays\n\n";
print FILE_SAVEOLD "  implicit none\n\n";
print FILE_SAVEOLD "  logical contains\n\n";
print FILE_SAVEOLD "  integer l,i\n\n";

# Write beginning of file simpleboundary.f90

print FILE_SIMPLEBOUNDARY "! Automatically generated file.  Do not edit!\n\n";
print FILE_SIMPLEBOUNDARY "  subroutine simpleboundary(l)\n\n";
print FILE_SIMPLEBOUNDARY "  use param\n";
print FILE_SIMPLEBOUNDARY "  use arrays\n\n";
print FILE_SIMPLEBOUNDARY "  implicit none\n\n";
print FILE_SIMPLEBOUNDARY "  logical contains\n\n";
print FILE_SIMPLEBOUNDARY "  integer l\n\n";

# Write beginning of file symmetries.f90

print FILE_SYMMETRIES "! Automatically generated file.  Do not edit!\n\n";
print FILE_SYMMETRIES "  subroutine symmetries(l)\n\n";
print FILE_SYMMETRIES "  use param\n";
print FILE_SYMMETRIES "  use arrays\n\n";
print FILE_SYMMETRIES "  implicit none\n\n";
print FILE_SYMMETRIES "  logical contains\n\n";
print FILE_SYMMETRIES "  integer i,l\n\n";
print FILE_SYMMETRIES "  do i=1,ghost\n\n";

# Write beginning of file syncall.f90

print FILE_SYNCALL "! Automatically generated file.  Do not edit!\n\n";
print FILE_SYNCALL "  subroutine syncall(l)\n\n";
print FILE_SYNCALL "  use param\n";
print FILE_SYNCALL "  use arrays\n\n";
print FILE_SYNCALL "  implicit none\n\n";
print FILE_SYNCALL "  logical contains\n\n";
print FILE_SYNCALL "  integer l\n\n";

# Write beginning of file update.f90

print FILE_UPDATE "! Automatically generated file.  Do not edit!\n\n";
print FILE_UPDATE "  subroutine update(l,dtw)\n\n";
print FILE_UPDATE "  use param\n";
print FILE_UPDATE "  use arrays\n\n";
print FILE_UPDATE "  implicit none\n\n";
print FILE_UPDATE "  logical contains\n\n";
print FILE_UPDATE "  integer l\n\n";
print FILE_UPDATE "  real(8) dtw \n\n";

# Write beginning of file boundinterp.inc

print FILE_BOUNDINTERP "! Automatically generated file.  Do not edit!\n\n";

# Write beginning of file restrict_copy.inc

print FILE_RESTRICTCOPY "! Automatically generated file.  Do not edit!\n\n";

# Write beginning of file restrict_send.inc

print FILE_RESTRICTSEND "! Automatically generated file.  Do not edit!\n\n";

# Write beginning of file restrict_recv.inc

print FILE_RESTRICTRECV "! Automatically generated file.  Do not edit!\n\n";

# Write beginning of file restrict_cosmo.inc

print FILE_RESTRICTCOSMO "! Automatically generated file.  Do not edit!\n\n";

# Open infile arrays.f90

open(INFILE,"src/base/arrays.config") or die "Can't open arrays.config: $!";

# Parse file arrays.f90 to identify declared arrays

my $line = " ";
my $nline = 0;

while ($line=<INFILE>) {

   $nline = $nline+1;

#  Look only for lines declaring real or complex arrays, ignore all other lines.

   if (($line =~ /^\s*REAL/i)||($line =~ /^\s*COMPLEX/i)) {

#     Check that all keywords are present and grab array name (make sure to
#     ignore possible comment at the end).

      if ($line =~ /^\s+REAL/i) {
          $type = "real(8)";
          $zero = "0.d0";
          if ($line =~ /REAL\s*(\S+)\s*!\s*SYMMETRY\s*=\s*(\S+)\s*,\s*INTENT\s*=\s*(\S+)\s*,\s*STORAGE\s*=\s*(.+)/i) {
             $var = $1;
             $sym = $2;
             $intent = $3;
             $storage = $4;
          } else {
             die "arrays.pl: Bad syntax for REAL array assignment in line ",$nline," of file arrays.config\n\n";
          }
      } elsif ($line =~ /^\s+COMPLEX/i) {
          $type = "complex(8)";
          $zero = "(0.d0,0.d0)";
          if ($line =~ /COMPLEX\s*(\S+)\s*!\s*SYMMETRY\s*=\s*(\S+)\s*,\s*INTENT\s*=\s*(\S+)\s*,\s*STORAGE\s*=\s*(.+)/i) {
	     $var = $1;
             $sym = $2;
             $intent = $3;
             $storage = $4;
          } else {
             die "arrays.pl: Bad syntax for COMPLEX array assignment in line ",$nline," of file arrays.config\n\n";
          }
      }

#     Check if we have an array that lives on phase space.

      if ($line =~ /PHASE/i) {
	  $phase = "true";
      } else {
	  $phase = "false";
      }

#     Check if we have a 0D array.

      if ($line =~ /ZEROD/i) {
	  $zerod = "true";
      } else {
	  $zerod = "false";
      }

#     Check if we have a cosmological array.

      if ($line =~ /cosmic_run/i) {
	  $cosmo = "true";
      } else {
	  $cosmo = "false";
      }

#     Check if the array has only one grid level.

      if ($line =~ /ONELEVEL/i) {
	  $onelevel = "true";
      } else {
	  $onelevel = "false";
      }

#     Check if the array does not require boundary conditions.

      if ($line =~ /NOBOUND/i) {
	  $nobound = "true";
      } else {
	  $nobound = "false";
      }

#     Check if the array should be checkpointed.

      if ($line =~ /CHECKPOINT/i) {
	  $checkpoint = "true";
      } else {
	  $checkpoint = "false";
      }

#     Check for commas.

      if ($var =~ /,/) {
          die "arrays.pl: Bad array assignment in line ",$nline," of file arrays.config\n",
              "           Only one array can be declared per line\n\n";
      }

#     Now write to FILE_ARRAYS code to declare arrays.

      if ($zerod eq "false") {

         if ($intent =~ /^POINTER$/i) {
            if ($onelevel eq "true") {
               print FILE_ARRAYS  "  ",$type,", dimension (:), pointer :: ",$var,"\n";
            } else {
               print FILE_ARRAYS  "  ",$type,", dimension (:,:), pointer :: ",$var,"\n";
            }
         } else {
            if ($onelevel eq "true") {
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:), target :: ",$var,"\n";
            } elsif ($intent =~ /^EVOLVE$/i) {
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:,:), target :: ",$var,"\n";
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:,:), target :: ",$var,"_p\n";
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:,:), target :: ",$var,"_a\n";
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:,:), target :: s",$var,"\n";
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:,:,:), target :: ",$var,"_bound\n";
            } else {
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:,:), target :: ",$var,"\n";
            }
         }

      } elsif ($zerod eq "true") {

         if ($intent =~ /^POINTER$/i) {
            print FILE_ARRAYS  "  ",$type,", dimension (:), pointer :: ",$var,"\n";
         } else {
            print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:), target :: ",$var,"\n";
            if ($intent =~ /^EVOLVE$/i) {
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:), target :: ",$var,"_p\n";
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:), target :: ",$var,"_a\n";
               print FILE_ARRAYS  "  ",$type,", allocatable, dimension (:), target :: s",$var,"\n";
            }
         }

      }

#     Write to FILE_ACCUMULATE code to save old variables.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_ACCUMULATE  "  if (",$cond,") then\n";
               print FILE_ACCUMULATE  "     if (k==1) then\n";
               print FILE_ACCUMULATE  "        ",$var,"_a(l,:) = w*s",$var,"(l,:)\n";
               print FILE_ACCUMULATE  "     else if (k<niter) then\n";
               print FILE_ACCUMULATE  "        ",$var,"_a(l,:) = ",$var,"_a(l,:) + w*s",$var,"(l,:)\n";
               print FILE_ACCUMULATE  "     else\n";
               print FILE_ACCUMULATE  "        s",$var,"(l,:)  = ",$var,"_a(l,:) + w*s",$var,"(l,:)\n";
               print FILE_ACCUMULATE  "     end if\n";
               print FILE_ACCUMULATE  "  end if\n\n";
            } else {
               print FILE_ACCUMULATE  "  if (k==1) then\n";
               print FILE_ACCUMULATE  "     ",$var,"_a(l,:) = w*s",$var,"(l,:)\n";
               print FILE_ACCUMULATE  "  else if (k<niter) then\n";
               print FILE_ACCUMULATE  "     ",$var,"_a(l,:) = ",$var,"_a(l,:) + w*s",$var,"(l,:)\n";
               print FILE_ACCUMULATE  "  else\n";
               print FILE_ACCUMULATE  "     s",$var,"(l,:)  = ",$var,"_a(l,:) + w*s",$var,"(l,:)\n";
               print FILE_ACCUMULATE  "  end if\n\n";
            }
         }

      } elsif ($zerod eq "true") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_ACCUMULATE  "  if (",$cond,") then\n";
               print FILE_ACCUMULATE  "     if (k==1) then\n";
               print FILE_ACCUMULATE  "        ",$var,"_a(l) = w*s",$var,"(l)\n";
               print FILE_ACCUMULATE  "     else if (k<niter) then\n";
               print FILE_ACCUMULATE  "        ",$var,"_a(l) = ",$var,"_a(l) + w*s",$var,"(l)\n";
               print FILE_ACCUMULATE  "     else\n";
               print FILE_ACCUMULATE  "        s",$var,"(l)  = ",$var,"_a(l) + w*s",$var,"(l)\n";
               print FILE_ACCUMULATE  "     end if\n";
               print FILE_ACCUMULATE  "  end if\n\n";
            } else {
               print FILE_ACCUMULATE  "  if (k==1) then\n";
               print FILE_ACCUMULATE  "     ",$var,"_a(l) = w*s",$var,"(l)\n";
               print FILE_ACCUMULATE  "  else if (k<niter) then\n";
               print FILE_ACCUMULATE  "     ",$var,"_a(l) = ",$var,"_a(l) + w*s",$var,"(l)\n";
               print FILE_ACCUMULATE  "  else\n";
               print FILE_ACCUMULATE  "     s",$var,"(l)  = ",$var,"_a(l) + w*s",$var,"(l)\n";
               print FILE_ACCUMULATE  "  end if\n\n";
            }
         }
      }

#     Write to FILE_ALLOCATEARRAYS code to allocate memory.

      if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
          $cond = $1;
          $space   = "   ";
          $newline = "";
          print FILE_ALLOCATEARRAYS  "  if (",$cond,") then\n";
      } elsif ($storage !~ /^ALWAYS/i) {
          print $storage,"\n";
          die "arrays.pl: Bad STORAGE assignment in line ",$nline," of file arrays.config\n\n";
      } else {
          $space   = "";
          $newline = "\n";
      }

      if ($zerod eq "false") {

         if ($intent =~/^OUTPUT$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (contains(outvars0D,\"",$var,"\").or. &\n";
            print FILE_ALLOCATEARRAYS  $space,"      contains(outvars1D,\"",$var,"\")) then\n";
            print FILE_ALLOCATEARRAYS  $space,"     if (trim(status)=='on') then\n";
            print FILE_ALLOCATEARRAYS  $space,"        allocate(",$var,"(0:Nl-1,1-ghost:Nrmax))\n";
            print FILE_ALLOCATEARRAYS  $space,"        ",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     else\n";
            print FILE_ALLOCATEARRAYS  $space,"        deallocate(",$var,")\n";
            print FILE_ALLOCATEARRAYS  $space,"     end if\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^EVOLVE$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (trim(status)=='on') then\n";
            print FILE_ALLOCATEARRAYS  $space,"     checkvars = trim(checkvars) // ',",$var,"'\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"(0:Nl-1,1-ghost:Nrmax)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(s",$var,"(0:Nl-1,1-ghost:Nrmax)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"_p(0:Nl-1,1-ghost:Nrmax)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"_a(0:Nl-1,1-ghost:Nrmax)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"_bound(0:Nl-1,0:ghost-1,0:3)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var,"   = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     s",$var,"  = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var,"_p = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var,"_a = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var,"_bound = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"  else\n";
            print FILE_ALLOCATEARRAYS  $space,"     deallocate(",$var,",s",$var,",",$var,"_p,",$var,"_a",",",$var,"_bound",")\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^AUXILIARY$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (trim(status)=='on') then\n";
            if ($checkpoint eq "true") {
               print FILE_ALLOCATEARRAYS  $space,"     checkvars = trim(checkvars) // ',",$var,"'\n";
            }
            if ($onelevel eq "true") {
               print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"(1-ghost:Nrmax))\n";
            } else {
               print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"(0:Nl-1,1-ghost:Nrmax))\n";
            }
            print FILE_ALLOCATEARRAYS  $space,"     ",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"  else\n";
            print FILE_ALLOCATEARRAYS  $space,"     deallocate(",$var,")\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^POINTER$/i) {
         } else {
            die "arrays.pl: Bad INTENT assignment in line ",$nline," of file arrays.config\n\n";
         }

      } elsif ($zerod eq "true") {

         if ($intent =~/^OUTPUT$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (contains(outvars0D,\"",$var,"\").or. &\n";
            print FILE_ALLOCATEARRAYS  $space,"      contains(outvars1D,\"",$var,"\")) then\n";
            print FILE_ALLOCATEARRAYS  $space,"     if (trim(status)=='on') then\n";
            print FILE_ALLOCATEARRAYS  $space,"        allocate(",$var,"(0:Nl-1))\n";
            print FILE_ALLOCATEARRAYS  $space,"        ",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     else\n";
            print FILE_ALLOCATEARRAYS  $space,"        deallocate(",$var,")\n";
            print FILE_ALLOCATEARRAYS  $space,"     end if\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^EVOLVE$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (trim(status)=='on') then\n";
            print FILE_ALLOCATEARRAYS  $space,"     checkvars = trim(checkvars) // ',",$var,"'\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"(0:Nl-1)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(s",$var,"(0:Nl-1)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"_p(0:Nl-1)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"_a(0:Nl-1)",")\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var,"   = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     s",$var,"  = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var,"_p = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var,"_a = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"  else\n";
            print FILE_ALLOCATEARRAYS  $space,"     deallocate(",$var,",s",$var,",",$var,"_p,",$var,"_a",")\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^AUXILIARY$/i) {
            print FILE_ALLOCATEARRAYS  $space,"  if (trim(status)=='on') then\n";
            if ($checkpoint eq "true") {
               print FILE_ALLOCATEARRAYS  $space,"     checkvars = trim(checkvars) // ',",$var,"'\n";
            }
            print FILE_ALLOCATEARRAYS  $space,"     allocate(",$var,"(0:Nl-1))\n";
            print FILE_ALLOCATEARRAYS  $space,"     ",$var," = ",$zero,"\n";
            print FILE_ALLOCATEARRAYS  $space,"  else\n";
            print FILE_ALLOCATEARRAYS  $space,"     deallocate(",$var,")\n";
            print FILE_ALLOCATEARRAYS  $space,"  end if\n",$newline;
         } elsif ($intent =~ /^POINTER$/i) {
         } else {
            die "arrays.pl: Bad INTENT assignment in line ",$nline," of file arrays.config\n\n";
         }

      }

      if ($storage =~ /^CONDITIONAL\s*\(.*\)/i) {
         print FILE_ALLOCATEARRAYS  "  else if (contains(outvars0D,\"",$var,"\").or. &\n";
         print FILE_ALLOCATEARRAYS  "           contains(outvars1D,\"",$var,"\")) then\n";
         print FILE_ALLOCATEARRAYS  "     if (rank==0) then\n";
         print FILE_ALLOCATEARRAYS  "        print *\n";
         print FILE_ALLOCATEARRAYS  "        print *, 'Error in parfile: array ",$var," has no storage,'\n";
         print FILE_ALLOCATEARRAYS  "        print *, 'so no output for it is possible.'\n";
         print FILE_ALLOCATEARRAYS  "        print *, 'Aborting! (subroutine allocatearrays.f90)'\n";
         print FILE_ALLOCATEARRAYS  "        print *\n";
         print FILE_ALLOCATEARRAYS  "     end if\n";
         print FILE_ALLOCATEARRAYS  "     call die\n";
         print FILE_ALLOCATEARRAYS  "  end if\n\n";
      }

#     Now write to FILE_GRABARRAY code to compare array name.

      if ($zerod eq "false") {

         if ($intent !~ /^POINTER$/i) {
            print FILE_GRABARRAY "  if (varname=='",$var,"') then\n";
            print FILE_GRABARRAY "     exists = .true.\n";
            print FILE_GRABARRAY "     savevar => ",$var,"\n";
            if ($intent =~ /^EVOLVE$/i) {
               print FILE_GRABARRAY "     grabvar_bound => ",$var,"_bound\n";
            }
            print FILE_GRABARRAY "  end if\n\n";
         }

         if ($intent =~ /^EVOLVE$/i) {
            print FILE_GRABARRAY "  if (varname=='s",$var,"') then\n";
            print FILE_GRABARRAY "     exists = .true.\n";
            print FILE_GRABARRAY "     savevar => s",$var,"\n";
            print FILE_GRABARRAY "  end if\n\n";
         }

      } elsif ($zerod eq "true") {

         if ($intent !~ /^POINTER$/i) {
            print FILE_GRABARRAY "  if (varname=='",$var,"') then\n";
            print FILE_GRABARRAY "     exists = .true.\n";
            print FILE_GRABARRAY "     savevar0D => ",$var,"\n";
            print FILE_GRABARRAY "     nullify(savevar)\n";
            print FILE_GRABARRAY "  end if\n\n";
         }

      }

#     Write to FILE_SAVEOLD code to save old variables.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_SAVEOLD  "  if (",$cond,") then\n";
               print FILE_SAVEOLD  "     ",$var,"_p(l,:) = ",$var,"(l,:)\n";
               print FILE_SAVEOLD  "     do i=0,ghost-1\n";
               print FILE_SAVEOLD  "        ",$var,"_bound(l,i,3) = ",$var,"_bound(l,i,2)\n";
               print FILE_SAVEOLD  "        ",$var,"_bound(l,i,2) = ",$var,"_bound(l,i,1)\n";
               print FILE_SAVEOLD  "        ",$var,"_bound(l,i,1) = ",$var,"(l,Nr-i)\n";
               print FILE_SAVEOLD  "     end do\n";
               print FILE_SAVEOLD  "  end if\n\n";
            } else {
               print FILE_SAVEOLD  "  ",$var,"_p(l,:) = ",$var,"(l,:)\n";
               print FILE_SAVEOLD  "  do i=0,ghost-1\n";
               print FILE_SAVEOLD  "     ",$var,"_bound(l,i,3) = ",$var,"_bound(l,i,2)\n";
               print FILE_SAVEOLD  "     ",$var,"_bound(l,i,2) = ",$var,"_bound(l,i,1)\n";
               print FILE_SAVEOLD  "     ",$var,"_bound(l,i,1) = ",$var,"(l,Nr-i)\n";
               print FILE_SAVEOLD  "  end do\n\n";
            }
         }

      } elsif ($zerod eq "true") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_SAVEOLD  "  if (",$cond,") then\n";
               print FILE_SAVEOLD  "     ",$var,"_p(l) = ",$var,"(l)\n";
               print FILE_SAVEOLD  "  end if\n\n";
            } else {
               print FILE_SAVEOLD  "  ",$var,"_p(l) = ",$var,"(l)\n";
            }
         }

      }

#     Write to FILE_SIMPLEBOUNDARIES code to apply simple boundary conditions.

      if ($zerod eq "false") {

         if (($intent =~ /EVOLVE/i) && ($nobound eq "false")) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
                $cond = $1;
                print FILE_SIMPLEBOUNDARY  "  if (",$cond,") then\n";
	        print FILE_SIMPLEBOUNDARY  "     if (boundtype=='static') then\n";
	        print FILE_SIMPLEBOUNDARY  "        s",$var,"(l,Nr) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "     else if (boundtype=='flat') then\n";
                print FILE_SIMPLEBOUNDARY  "        s",$var,"(l,Nr) = s",$var,"(l,Nr-1)\n";
                print FILE_SIMPLEBOUNDARY  "     end if\n";
                print FILE_SIMPLEBOUNDARY  "  end if\n\n";
            } else {
	        print FILE_SIMPLEBOUNDARY  "  if (boundtype=='static') then\n";
	        print FILE_SIMPLEBOUNDARY  "     s",$var,"(l,Nr) = 0.d0\n";
	        print FILE_SIMPLEBOUNDARY  "  else if (boundtype=='flat') then\n";
	        print FILE_SIMPLEBOUNDARY  "     s",$var,"(l,Nr) = s",$var,"(l,Nr-1)\n";
                print FILE_SIMPLEBOUNDARY  "  end if\n\n";
            }
         }

      }

#     Write to FILE_SYMMETRIES code to apply symmetries at origin.
#     Here I allow the possibility that the symmetry is not fixed
#     to 1 or -1, but is given by an expression that should evaluate
#     to +-1 depending on certain parameters.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_SYMMETRIES  "     if (",$cond,") then\n";
               if ($sym == "+1") {
                  print FILE_SYMMETRIES  "        ",$var,"(l,1-i) = + ",$var,"(l,i)\n";
               } elsif ($sym == "-1") {
                  print FILE_SYMMETRIES  "        ",$var,"(l,1-i) = - ",$var,"(l,i)\n";
               } else {
                  print FILE_SYMMETRIES  "        ",$var,"(l,1-i) = ",$sym,"*",$var,"(l,i)\n";
               }
               print FILE_SYMMETRIES  "     end if\n\n";
            } else {
               if ($sym == "+1") {
	          print FILE_SYMMETRIES  "     ",$var,"(l,1-i) = + ",$var,"(l,i)\n\n";
               } elsif ($sym == "-1") {
                  print FILE_SYMMETRIES  "     ",$var,"(l,1-i) = - ",$var,"(l,i)\n\n";
               } else {
                  print FILE_SYMMETRIES  "     ",$var,"(l,1-i) = ",$sym,"*",$var,"(l,i)\n\n";
               }
            }
         }

      }

#     Write to FILE_SYNCALL code to synchronize across processors.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
                $cond = $1;
                print FILE_SYNCALL  "  if (",$cond,") then\n";
	        print FILE_SYNCALL  "     syncvar => ",$var,"(l,:)\n";
	        print FILE_SYNCALL  "     call sync\n";
                print FILE_SYNCALL  "  end if\n\n";
            } else {
	        print FILE_SYNCALL  "  syncvar => ",$var,"(l,:)\n";
	        print FILE_SYNCALL  "  call sync\n\n";
            }
         }

      }

#     Write to FILE_UPDATE code to update variables.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_UPDATE  "  if (",$cond,") then\n";
               print FILE_UPDATE  "     ",$var,"(l,:) = ",$var,"_p(l,:) + dtw*s",$var,"(l,:)\n";
               print FILE_UPDATE  "  end if\n\n";
            } else {
               print FILE_UPDATE  "  ",$var,"(l,:) = ",$var,"_p(l,:) + dtw*s",$var,"(l,:)\n\n";
            }
         }

      } elsif ($zerod eq "true") {

        if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_UPDATE  "  if (",$cond,") then\n";
               print FILE_UPDATE  "     ",$var,"(l) = ",$var,"_p(l) + dtw*s",$var,"(l)\n";
               print FILE_UPDATE  "  end if\n\n";
            } else {
               print FILE_UPDATE  "  ",$var,"(l) = ",$var,"_p(l) + dtw*s",$var,"(l)\n\n";
            }
        }

      }

#     Write to FILE_BOUNDINTERP code to interpolate variables at boundaries.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
         #if (($intent =~ /EVOLVE/i) && ($nobound eq "false")) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_BOUNDINTERP  "  if (",$cond,") then\n";
               print FILE_BOUNDINTERP  "     interpvar => ",$var,"\n";
               print FILE_BOUNDINTERP  "     aux1 = interp(l-1,r0,.false.)\n";
               print FILE_BOUNDINTERP  "     call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
               print FILE_BOUNDINTERP  "     if (border==1) then\n";
               print FILE_BOUNDINTERP  "        ",$var,"(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*",$var,"_p(l,Nr-i)\n";
               print FILE_BOUNDINTERP  "     else\n";
               print FILE_BOUNDINTERP  "        ",$var,"(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
               print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound(l,i,1) &\n";
               print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound(l,i,2)\n";
               print FILE_BOUNDINTERP  "     end if\n";
               print FILE_BOUNDINTERP  "  end if\n\n";
            } else {
               print FILE_BOUNDINTERP  "  interpvar => ",$var,"\n";
               print FILE_BOUNDINTERP  "  aux1 = interp(l-1,r0,.false.)\n";
               print FILE_BOUNDINTERP  "  call MPI_ALLREDUCE(aux1,aux2,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)\n";
               print FILE_BOUNDINTERP  "  if (border==1) then\n";
               print FILE_BOUNDINTERP  "     ",$var,"(l,Nr-i) = (tl-t0)/(tp-t0)*aux2 + (tl-tp)/(t0-tp)*",$var,"_p(l,Nr-i)\n";
               print FILE_BOUNDINTERP  "  else\n";
               print FILE_BOUNDINTERP  "     ",$var,"(l,Nr-i) = (tl-t0)*(tl-tm1)/((tp-t0)*(tp-tm1))*aux2 &\n";
               print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-tm1)/((t0 -tp)*(t0-tm1))*",$var,"_bound(l,i,1) &\n";
               print FILE_BOUNDINTERP  "              + (tl-tp)*(tl-t0 )/((tm1-tp)*(tm1-t0))*",$var,"_bound(l,i,2)\n";
               print FILE_BOUNDINTERP  "  end if\n\n";
            }
         }

      }

#     Write to FILE_RESTRICTCOPY code to restrict on processor 0.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               $cond = $1;
               print FILE_RESTRICTCOPY  "  if (",$cond,") then\n";
               print FILE_RESTRICTCOPY  "     do i=1,Nr-ghost,2\n";
               print FILE_RESTRICTCOPY  "        r0 = r(l,i) + delta\n";
               print FILE_RESTRICTCOPY  "        interpvar => ",$var,"\n";
               print FILE_RESTRICTCOPY  "        ",$var,"(l-1,i/2+1) = interp(l,r0,.true.)\n";
               print FILE_RESTRICTCOPY  "     end do\n";
               print FILE_RESTRICTCOPY  "  end if\n\n";
            } else {
               print FILE_RESTRICTCOPY  "  do i=1,Nr-ghost,2\n";
               print FILE_RESTRICTCOPY  "     r0 = r(l,i) + delta\n";
               print FILE_RESTRICTCOPY  "     interpvar => ",$var,"\n";
               print FILE_RESTRICTCOPY  "     ",$var,"(l-1,i/2+1) = interp(l,r0,.true.)\n";
               print FILE_RESTRICTCOPY  "  end do\n\n";
           }
         }

      }

#     Write to FILE_RESTRICTSEND code to send data.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
           if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               print FILE_RESTRICTSEND  "  if (",$cond,") then\n";
               print FILE_RESTRICTSEND  "     do i=imin,Nr-ghost,2\n";
               print FILE_RESTRICTSEND  "        r0 = r(l,i) + delta\n";
               print FILE_RESTRICTSEND  "        interpvar => ",$var,"\n";
               print FILE_RESTRICTSEND  "        w(i/2) = interp(l,r0,.true.)\n";
               print FILE_RESTRICTSEND  "     end do\n";
               print FILE_RESTRICTSEND  "     call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)\n";
               print FILE_RESTRICTSEND  "  end if\n\n";
            } else {
               print FILE_RESTRICTSEND  "  do i=imin,Nr-ghost,2\n";
               print FILE_RESTRICTSEND  "     r0 = r(l,i) + delta\n";
               print FILE_RESTRICTSEND  "     interpvar => ",$var,"\n";
               print FILE_RESTRICTSEND  "     w(i/2) = interp(l,r0,.true.)\n";
               print FILE_RESTRICTSEND  "  end do\n";
               print FILE_RESTRICTSEND  "  call MPI_SEND(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,ierr)\n\n";
            }
         }

      }

#     Write to FILE_RESTRICTRECV code to receive data.

      if ($zerod eq "false") {

         if ($intent =~ /EVOLVE/i) {
            if ($storage =~ /^CONDITIONAL\s*\((.*)\)/i) {
               print FILE_RESTRICTRECV  "  if (",$cond,") then\n";
               print FILE_RESTRICTRECV  "     call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)\n";
               print FILE_RESTRICTRECV  "     do i=imin,Nrl(p)-ghost,2\n";
               print FILE_RESTRICTRECV  "        ",$var,"(l-1,i/2+k) = w(i/2)\n";
               print FILE_RESTRICTRECV  "     end do\n";
               print FILE_RESTRICTRECV  "  end if\n\n";
            } else {
               print FILE_RESTRICTRECV  "  call MPI_RECV(w,Ndata,MPI_REAL8,p,1,MPI_COMM_WORLD,status,ierr)\n";
               print FILE_RESTRICTRECV  "  do i=imin,Nrl(p)-ghost,2\n";
               print FILE_RESTRICTRECV  "     ",$var,"(l-1,i/2+k) = w(i/2)\n";
               print FILE_RESTRICTRECV  "  end do\n\n";
            }
         }

      }

#     Write to FILE_RESTRICTCOSMO code to restrict background cosmological functions.

      if ($zerod eq "true"  && $cosmo eq "true") {
         if ($intent =~ /EVOLVE/i) {
            print FILE_RESTRICTCOSMO  "  ",$var,"(l-1) = ",$var,"(l)\n\n";
         }
      }

#  Close two main conditional statements.

   } elsif (($line !~ /^\s*!.*/)&&($line !~ /^\s*$/)) { 
      die "arrays.pl: Bad syntax in line ",$nline," of file arrays.config\n\n";
   }
}

# Close INFILE.

close(INFILE);

# Write ending of file arrays.f90.

print FILE_ARRAYS "\n  end module arrays\n\n";

# Write ending of file accumulate.f90.

print FILE_ACCUMULATE "  end subroutine accumulate\n\n";

# Write ending of file allocatearrays.f90

print FILE_ALLOCATEARRAYS "  end subroutine allocatearrays\n\n";

# Write ending of file grabarray.f90.

print FILE_GRABARRAY "  if (.not.exists) then\n";
print FILE_GRABARRAY "     if (rank==0) then\n";
print FILE_GRABARRAY "        print *\n";
print FILE_GRABARRAY "        print *, 'Error in parfile, non-existent array: ',varname\n";
print FILE_GRABARRAY "        print *, 'Aborting! (subroutine grabarray.f90)'\n";
print FILE_GRABARRAY "        print *\n";
print FILE_GRABARRAY "     end if\n";
print FILE_GRABARRAY "     call die\n";
print FILE_GRABARRAY "  end if\n\n";
print FILE_GRABARRAY "  end subroutine grabarray\n\n";

# Write ending of file saveold.f90.

print FILE_SAVEOLD "  end subroutine saveold\n\n";

# Write beginning of file simpleboundary.f90

print FILE_SIMPLEBOUNDARY "  end subroutine simpleboundary\n\n";

# Write ending of file symmetries.f90.

print FILE_SYMMETRIES  "  end do\n\n";
print FILE_SYMMETRIES  "  end subroutine symmetries\n\n";

# Write ending of file syncall.f90.

print FILE_SYNCALL  "  end subroutine syncall\n\n";

# Write ending of file update.f90.

print FILE_UPDATE "  end subroutine update\n\n";

# Close output files.

close(FILE_ARRAYS);
close(FILE_ACCUMULATE);
close(FILE_ALLOCATEARRAYS);
close(FILE_GRABARRAY);
close(FILE_SAVEOLD);
close(FILE_SIMPLEBOUNDARY);
close(FILE_SYMMETRIES);
close(FILE_SYNCALL);
close(FILE_UPDATE);

close(FILE_BOUNDINTERP);
close(FILE_RESTRICTCOPY);
close(FILE_RESTRICTSEND);
close(FILE_RESTRICTRECV);
close(FILE_RESTRICTCOSMO);

