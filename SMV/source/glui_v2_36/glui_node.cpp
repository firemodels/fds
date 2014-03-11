/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_node.cpp - linked-list tree structure


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  WWW:    http://sourceforge.net/projects/glui/
  Forums: http://sourceforge.net/forum/?group_id=92496

  This software is provided 'as-is', without any express or implied 
  warranty. In no event will the authors be held liable for any damages 
  arising from the use of this software. 

  Permission is granted to anyone to use this software for any purpose, 
  including commercial applications, and to alter it and redistribute it 
  freely, subject to the following restrictions: 

  1. The origin of this software must not be misrepresented; you must not 
  claim that you wrote the original software. If you use this software 
  in a product, an acknowledgment in the product documentation would be 
  appreciated but is not required. 
  2. Altered source versions must be plainly marked as such, and must not be 
  misrepresented as being the original software. 
  3. This notice may not be removed or altered from any source distribution. 

*****************************************************************************/

#include "GL/glui.h"
#include "glui_internal.h"

/********************************************* GLUI_Node::GLUI_Node() *******/

GLUI_Node::GLUI_Node()
: 
    parent_node(NULL),
    child_head(NULL),
    child_tail(NULL),
    next_sibling(NULL),
    prev_sibling(NULL)
{
}

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

/**************************************** GLUI_Node::dump() **************/

void GLUI_Node::dump( FILE *out, const char *name )
{
    fprintf( out, "GLUI_node: %s\n", name );
    fprintf( out, "   parent: %p     child_head: %p    child_tail: %p\n",
        (void *) parent_node,
        (void *) child_head,
        (void *) child_tail );
    fprintf( out, "   next: %p       prev: %p\n",
        (void *) next_sibling,
        (void *) prev_sibling );
}
