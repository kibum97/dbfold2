#ifndef CONTACTS_H
#define CONTACTS_H

#include "globals.h" /* for global variables */

void Contacts(struct Context *ctx, const struct Topology *top, const struct System *sys);
void TypeContacts(struct Context *ctx, const struct System *sys, const struct Topology *top);
void CheckForContacts(struct Context *ctx, const struct System *sys, short a, short b);
unsigned char CheckForDeltaContacts(struct contact_data *Data, struct int_vector XX,
                                    struct int_vector YY, short type_a, short type_b,
                                    const struct System *sys, struct monte_carlo_flags mc_flags);
void          check_bb_contacts(struct Context *ctx, short a, short b, const struct Topology *top,
                                const struct Simulation *sim);
void          NewDeltaContacts(struct Context *ctx, const struct System *sys, short rotate_natoms,
                               short *rotate_atom, char *not_rotated);

/* wmj */
void   fill_calpha_contact_string(const struct Topology *top, const struct residue *residues,
                                  const struct atom *atoms, short *contactstring);
double number_of_calpha_contacts(const struct Topology *top, const struct residue *residues,
                                 const struct atom *atoms);
double number_of_calpha_native_contacts(const struct Topology *top, const struct residue *residues,
                                        const struct atom *atoms, const short *contactstring);
double number_of_calpha_nonnative_contacts(const struct Topology *top,
                                           const struct residue *residues, const struct atom *atoms,
                                           const short *contactstring);
double hamming_distance_calpha_contacts(const struct Topology *top, const struct residue *residues,
                                        const struct atom *atoms, const short *contactstring);
/* end wmj */

#endif
