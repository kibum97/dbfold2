#ifndef CONTACTS_H
#define CONTACTS_H

#include "globals.h" /* for global variables */

void          Contacts();
void          TypeContacts();
void          DTypeContacts();
short         Bin(long);
void          CheckForContacts(short, short);
unsigned char CheckForDeltaContacts(struct contact_data *, struct int_vector, struct int_vector,
                                    short, short);
void          check_bb_contacts(short a, short b);
void          NewDeltaContacts(short rotate_natoms, short *rotate_atom, char *not_rotated);

/* wmj */
void fill_calpha_contact_string(const struct residue *residues, const struct atom *atoms,
                                  short *contactstring);
double number_of_calpha_contacts(const struct residue *residues, const struct atom *atoms);
double number_of_calpha_native_contacts(const struct residue *residues, const struct atom *atoms,
                                        const short *contactstring);
double number_of_calpha_nonnative_contacts(const struct residue *residues, const struct atom *atoms,
                                           const short *contactstring);
double hamming_distance_calpha_contacts(const struct residue *residues, const struct atom *atoms,
                                        const short *contactstring);
/* end wmj */

#endif

