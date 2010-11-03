/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_node.cpp - linked-list tree structure


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/


#include "glui.h"
#include "stdinc.h"

/********************************************* GLUI_Node::first() *******/
/* Returns first sibling in 'this' node's sibling list                  */

GLUI_Node   *GLUI_Node::first_sibling( void )
{
  if ( parent_node == NULL )  
    return this;           /* root node has no siblings */
  else
    return parent_node->child_head;
}


/******************************************** GLUI_Node::next() ********/
/* Returns next sibling in 'this' node's sibling list                  */

GLUI_Node    *GLUI_Node::next( void )
{
  return next_sibling;
}


/******************************************** GLUI_Node::prev() ********/
/* Returns prev sibling in 'this' node's sibling list                  */

GLUI_Node    *GLUI_Node::prev( void )
{
  return prev_sibling;
}


/********************************************* GLUI_Node::last() *******/
/* Returns last sibling in 'this' node's sibling list                  */

GLUI_Node   *GLUI_Node::last_sibling( void )
{
  if ( parent_node == NULL )
    return this;            /* root node has no siblings */
  else
    return parent_node->child_tail;
}


/*************************** GLUI_Node::link_this_to_parent_last() *******/
/* Links as last child of parent                                         */

void   GLUI_Node::link_this_to_parent_last( GLUI_Node *new_parent )
{
  if ( new_parent->child_tail == NULL ) {   /* parent has no children */
    new_parent->child_head = this;
    new_parent->child_tail = this;
    this->parent_node      = new_parent;
    this->parent_node      = new_parent;
  }
  else {                                 /* parent has children */
    new_parent->child_tail->next_sibling = this;
    this->prev_sibling                   = new_parent->child_tail;
    new_parent->child_tail               = this;
    this->parent_node                    = new_parent;
  }
}


/*************************** GLUI_Node::link_this_to_parent_first() *******/
/* Links as first child of parent                                         */

void   GLUI_Node::link_this_to_parent_first( GLUI_Node *new_parent )
{
  if ( new_parent->child_head == NULL ) {   /* parent has no children */
    new_parent->child_head               = this;
    new_parent->child_tail               = this;
    this->parent_node                    = new_parent;
  }
  else {                                 /* parent has children */
    new_parent->child_head->prev_sibling = this;
    this->next_sibling                   = new_parent->child_head;
    new_parent->child_head               = this;
    this->parent_node                    = new_parent;
  }
}

/**************************** GLUI_Node::link_this_to_sibling_next() *****/

void   GLUI_Node::link_this_to_sibling_next( GLUI_Node *sibling )
{
  if ( sibling->next_sibling == NULL ) {    /* node has no next sibling */
    sibling->next_sibling  = this;
    this->prev_sibling     = sibling;

    /* This was the parent's last child, so update that as well */
    if ( sibling->parent_node  != NULL ) {
      sibling->parent_node->child_tail = this;
    }
  }
  else {                            /* node already has a next sibling */
    sibling->next_sibling->prev_sibling = this;
    this->next_sibling                  = sibling->next_sibling;
    sibling->next_sibling               = this;
    this->prev_sibling                  = sibling;
  }

  this->parent_node = sibling->parent_node;
}


/**************************** GLUI_Node::link_this_to_sibling_prev() *****/

void   GLUI_Node::link_this_to_sibling_prev( GLUI_Node *sibling )
{
  if ( sibling->prev_sibling == NULL ) {    /* node has no prev sibling */
    sibling->prev_sibling  = this;
    this->next_sibling     = sibling;

    /* This was the parent's first child, so update that as well */
    if ( sibling->parent_node  != NULL ) {
      sibling->parent_node->child_head = this;
    }
  }
  else {                            /* node already has a prev sibling */
    sibling->prev_sibling->next_sibling = this;
    this->prev_sibling                  = sibling->prev_sibling;
    sibling->prev_sibling               = this;
    this->next_sibling                  = sibling;
  }

  this->parent_node = sibling->parent_node;
}

/**************************************** GLUI_Node::unlink() **************/

void   GLUI_Node::unlink( void )
{
  /* Unlink from prev sibling */
  if ( this->prev_sibling != NULL ) {
    this->prev_sibling->next_sibling = this->next_sibling;
  }
  else {                 /* No prev sibling: this was parent's first child */
    this->parent_node->child_head = this->next_sibling;
  }

  /* Unlink from next sibling */
  if ( this->next_sibling != NULL ) {
    this->next_sibling->prev_sibling = this->prev_sibling;
  }
  else {                /* No next sibling: this was parent's last child */
    this->parent_node->child_tail = this->prev_sibling;
  }

  this->parent_node  = NULL;
  this->next_sibling = NULL;
  this->prev_sibling = NULL;
  this->child_head   = NULL;
  this->child_tail   = NULL;
}
