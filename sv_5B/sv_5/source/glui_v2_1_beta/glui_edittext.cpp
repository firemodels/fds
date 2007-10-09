/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_edittext.cpp - GLUI_EditText control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/****************************** GLUI_EditText::mouse_down_handler() **********/

int    GLUI_EditText::mouse_down_handler( int local_x, int local_y )
{
  int tmp_insertion_pt;

  if ( debug )    dump( stdout, "-> MOUSE DOWN" );

  tmp_insertion_pt = find_insertion_pt( local_x, local_y );  
  if ( tmp_insertion_pt == -1 ) {
    if ( glui )
      glui->disactivate_current_control(  );
    return false;
  }

  insertion_pt = tmp_insertion_pt;

  sel_start = sel_end = insertion_pt;

  if ( can_draw())
    update_and_draw_text();

  if ( debug )    dump( stdout, "<- MOUSE UP" );

  return true;
}


/******************************** GLUI_EditText::mouse_up_handler() **********/

int    GLUI_EditText::mouse_up_handler( int local_x, int local_y, int inside )
{
  return false;
}


/***************************** GLUI_EditText::mouse_held_down_handler() ******/

int    GLUI_EditText::mouse_held_down_handler( int local_x, int local_y,
					       int new_inside)
{
  int tmp_pt;

  if ( NOT new_inside ) 
    return false;

  if ( debug )    dump( stdout, "-> HELD DOWN" );
  
  tmp_pt = find_insertion_pt( local_x, local_y );
  
  if ( tmp_pt == -1 AND sel_end != 0 ) {    /* moved mouse past left edge */
    special_handler( GLUT_KEY_LEFT, GLUT_ACTIVE_SHIFT );
  }
  else if ( tmp_pt == substring_end+1 AND sel_end != (int) strlen(text)) {    
    /* moved mouse past right edge */
    special_handler( GLUT_KEY_RIGHT, GLUT_ACTIVE_SHIFT );    
  }
  else if ( tmp_pt != -1 AND tmp_pt != sel_end ) {
    sel_end = insertion_pt = tmp_pt;
    
    update_and_draw_text();
  }

  if ( debug )
    dump( stdout, "<- HELD DOWN" );

  return false;
}


/****************************** GLUI_EditText::key_handler() **********/

int    GLUI_EditText::key_handler( unsigned char key,int modifiers )
{
  int i, regular_key;
  /* int has_selection;              */

  if ( NOT glui )
    return false;

  if ( debug )
    dump( stdout, "-> KEY HANDLER" );

  regular_key = false;
  /*  has_selection = (sel_start != sel_end);              */

  if ( key == 21 AND (modifiers & GLUT_ACTIVE_CTRL )!=0) { /* DEL all text */
    /** This one (key==21) may not port!! */
    
    insertion_pt = -1;  
    text[0] = '\0';
    sel_start = sel_end = -1;
  }
  else if ( key == 13 ) {           /* RETURN */
    /*    glui->disactivate_current_control();              */
    disactivate();  /** Force callbacks, etc **/
    activate(GLUI_ACTIVATE_TAB);     /** Reselect all text **/
    translate_and_draw_front();
    return true;
  }
  else if ( key  == 27 ) {         /* ESCAPE */
    glui->disactivate_current_control();
    return true;
  }
  else if ( key == 8 ) {       /* BACKSPACE */
    if ( sel_start == sel_end ) {   /* no selection */
      if ( insertion_pt > 0 ) {
	/*** See if we're deleting a period in a float data-type box ***/
	if ( data_type == GLUI_EDITTEXT_FLOAT AND text[insertion_pt-1]=='.' )
	  num_periods--;
	
	/*** Shift over string first ***/
	insertion_pt--;
	for( i=insertion_pt; i< (int)strlen( text ); i++ )
	  text[i] = text[i+1];    
      }
    }
    else {                         /* There is a selection */
      clear_substring( MIN(sel_start,sel_end), MAX(sel_start,sel_end ));
      insertion_pt = MIN(sel_start,sel_end);
      sel_start = sel_end = insertion_pt;
    }
  }
  else {                      /* Regular key */    
    regular_key = true;
    
    /** Check if we only accept numbers **/
    if (data_type == GLUI_EDITTEXT_FLOAT ) {
      if ( (key < '0' OR key > '9') AND key != '.' AND key != '-' )
	return true;

      if ( key == '-' ) { /* User typed a '-' */

	/* If user has first character selected, then '-' is allowed */
	if ( NOT ( MIN(sel_start,sel_end) == 0 AND
		   MAX(sel_start,sel_end) > 0 ) ) {

	  /* User does not have 1st char selected */
	  if (insertion_pt != 0 OR text[0] == '-' ) {
	    return true; /* Can only place negative at beginning of text,
			    and only one of them */
	  }
	}
      }

      if ( key == '.' ) {
	/*printf( "PERIOD: %d\n", num_periods );              */

	if ( num_periods > 0 ) {
	  /** We're trying to type a period, but the text already contains
	    a period.  Check whether the period is contained within
	    is current selection (thus it will be safely replaced) **/

	  int period_found = false; 
	  if ( sel_start != sel_end ) {
	    for( i=MIN(sel_end,sel_start); i<MAX(sel_start,sel_end); i++ ) {
	      /*  printf( "%c ", text[i] );              */
	      if ( text[i] == '.' ) {
		period_found = true;
		break;
	      }
	    }
	  }

	  /* printf( "found: %d    num: %d\n", period_found, num_periods );              */
	  
	  if ( NOT period_found )
	    return true;
	}
      }
    } 
    else if (data_type == GLUI_EDITTEXT_INT)	
    {
      if ( (key < '0' OR key > '9') AND key != '-' )
	return true;

      if ( key == '-' ) { /* User typed a '-' */

	/* If user has first character selected, then '-' is allowed */
	if ( NOT ( MIN(sel_start,sel_end) == 0 AND
		   MAX(sel_start,sel_end) > 0 ) ) {

	  /* User does not have 1st char selected */
	  if (insertion_pt != 0 OR text[0] == '-' ) {
	    return true; /* Can only place negative at beginning of text,
			    and only one of them */
	  }
	}
      }
    }

    /** This is just to get rid of warnings - the flag regular_key is 
      set if the key was not a backspace, return, whatever.  But I
      believe if we're here, we know it was a regular key anyway */
    if ( regular_key ) {
    }

    /**** If there's a current selection, erase it ******/
    if ( sel_start != sel_end ) {
      clear_substring( MIN(sel_start,sel_end), MAX(sel_start,sel_end ));
      insertion_pt = MIN(sel_start,sel_end);
      sel_start = sel_end = insertion_pt;
    }

    /******** check whether we have space ******/
    if ( (int)strlen( text ) + 2 >= sizeof( GLUI_String ))
      return false;

    /******** We insert the character into the string ***/
     
    /*** Shift over string first ***/
    for( i=(int)strlen( text ); i >= insertion_pt; i-- )
      text[i+1] = text[i];
    
    /******** Now insert the character ********/
    text[insertion_pt] = key;    

    /******** Move the insertion point and substring_end one over ******/
    insertion_pt++;
    substring_end++;

    sel_start = sel_end = insertion_pt;
  }
  
  /******** Now redraw text ***********/
  /* Hack to prevent text box from being cleared first **/  
  /**  int substring_change =  update_substring_bounds();
    draw_text_only = 
    (NOT substring_change AND NOT has_selection AND regular_key ); 
    */

  draw_text_only = false;  /** Well, hack is not yet working **/
  update_and_draw_text();
  draw_text_only = false;


  if ( debug )
    dump( stdout, "<- KEY HANDLER" );

  /*** Now look to see if this string has a period ***/
  num_periods = 0;
  for( i=0; i<(int)strlen(text); i++ )
    if ( text[i] == '.' )
      num_periods++;

  return true;
}


/****************************** GLUI_EditText::activate() **********/

void    GLUI_EditText::activate( int how )
{
  if ( debug )
    dump( stdout, "-> ACTIVATE" );

  active = true;

  if ( how == GLUI_ACTIVATE_MOUSE )
    return;  /* Don't select everything if activated with mouse */

  strcpy( orig_text, text );

  sel_start    = 0;
  sel_end      = (int)strlen(text);
  insertion_pt = 0;

  if ( debug )
    dump( stdout, "<- ACTIVATE" );
}


/****************************** GLUI_EditText::disactivate() **********/

void    GLUI_EditText::disactivate( void )
{
  int    new_int_val;
  float  new_float_val;

  active = false;

  if ( NOT glui )
    return;

  if ( debug )
    dump( stdout, "-> DISACTIVATE" );

  sel_start = sel_end = insertion_pt = -1; 

  /***** Retrieve the current value from the text *****/
  /***** The live variable will be updated by set_text() ****/
  if ( data_type == GLUI_EDITTEXT_FLOAT ) {
    if ( text[0] == '\0' ) /* zero-length string - make it "0.0" */
      strcpy( text, "0.0" );

    new_float_val = atof( text );

    set_float_val( new_float_val );
  }
  else if ( data_type == GLUI_EDITTEXT_INT ) {
    if ( text[0] == '\0' ) /* zero-length string - make it "0" */
      strcpy( text, "0" );

    new_int_val = atoi( text );

    set_int_val( new_int_val );
  }
  else 
    if ( data_type == GLUI_EDITTEXT_TEXT ) {
      set_text(text); /* This will force callbacks and gfx refresh */
    }

  update_substring_bounds();

  /******** redraw text without insertion point ***********/
  translate_and_draw_front();

  /***** Now do callbacks if value changed ******/
  if ( strcmp( orig_text, text ) != 0 ) {
    this->execute_callback();
    
    if ( 0 ) {
      /* THE CODE BELOW IS FROM WHEN SPINNER ALSO MAINTAINED CALLBACKS    */
      if ( spinner == NULL ) {   /** Are we independent of a spinner?  **/  
	if ( callback ) {              
	  callback( this->user_id );              
	}              
      }              
      else {                      /* We're attached to a spinner */              
	spinner->do_callbacks();  /* Let the spinner do the callback stuff */  
      }              
    }
  }

  if ( debug )
    dump( stdout, "<- DISACTIVATE" );
}

/****************************** GLUI_EditText::draw() **********/

void    GLUI_EditText::draw( int x, int y )
{
  int orig;
  int name_x;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  name_x = MAX(text_x_offset - string_width(this->name) - 3,0);
  draw_name( name_x , 13);

  glBegin( GL_LINES );
  glColor3f( .5, .5, .5 );
  glVertex2i( text_x_offset, 0 );     glVertex2i( w, 0 );
  glVertex2i( text_x_offset, 0 );     glVertex2i( text_x_offset, h );     

  glColor3f( 1., 1., 1. );
  glVertex2i( text_x_offset, h );     glVertex2i( w, h );
  glVertex2i( w, h );                 glVertex2i( w, 0 );

  if ( enabled )
    glColor3f( 0., 0., 0. );
  else
    glColor3f( .25, .25, .25 );
  glVertex2i( text_x_offset+1, 1 );     glVertex2i( w-1, 1 );
  glVertex2i( text_x_offset+1, 1 );     glVertex2i( text_x_offset+1, h-1 );

  glColor3f( .75, .75, .75 );
  glVertex2i( text_x_offset+1, h-1 );     glVertex2i( w-1, h-1 );
  glVertex2i( w-1, h-1 );                 glVertex2i( w-1, 1 );
  glEnd();

  /** Find where to draw the text **/
  update_substring_bounds();
  draw_text(0,0);
  
  draw_insertion_pt();

  restore_window(orig);
}



/************************** GLUI_EditText::update_substring_bounds() *********/

int    GLUI_EditText::update_substring_bounds( void )
{
  int box_width;
  int text_len = (int)strlen(text);
  int old_start, old_end;

  old_start = substring_start;
  old_end = substring_end;

  /*** Calculate the width of the usable area of the edit box ***/
  box_width = MAX( this->w - this->text_x_offset 
		   - 4     /*  2 * the two-line box border */ 
		   - 2 * GLUI_EDITTEXT_BOXINNERMARGINX, 0 );

  CLAMP( substring_end, 0, MAX(text_len-1,0) );
  CLAMP( substring_start, 0, MAX(text_len-1,0) );

  if ( debug )    dump( stdout, "-> UPDATE SS" );

  if ( insertion_pt >= 0 AND 
       insertion_pt < substring_start ) {   /* cursor moved left */
    substring_start = insertion_pt;

    while ( substring_width( substring_start, substring_end ) > box_width )
      substring_end--;
  }
  else if ( insertion_pt > substring_end ) {  /* cursor moved right */
    substring_end = insertion_pt-1;

    while ( substring_width( substring_start, substring_end ) > box_width )
      substring_start++;
  }
  else {   /* cursor is within old substring bounds */
    if ( last_insertion_pt > insertion_pt ) {  /* cursor moved left */
    }
    else {
      while ( substring_width( substring_start, substring_end ) > box_width )
	substring_end--;

      while(substring_width( substring_start, substring_end+1 ) <= box_width
	    AND substring_end < text_len-1 )
      	substring_end++;
    }
  }

  while ( substring_width( substring_start, substring_end ) > box_width )
    substring_end--;

  last_insertion_pt = insertion_pt;

  /*** No selection if not enabled ***/
  if ( NOT enabled ) {
    sel_start = sel_end = 0;
  }

  if ( debug )    dump( stdout, "<- UPDATE SS" );

  if ( substring_start == old_start AND substring_end == old_end )
    return false;  /*** bounds did not change ***/
  else 
    return true;   /*** bounds did change ***/
}


/********************************* GLUI_EditText::update_x_offsets() *********/

void    GLUI_EditText::update_x_offsets( void )
{
}
 

/********************************* GLUI_EditText::draw_text() ****************/

void    GLUI_EditText::draw_text( int x, int y )
{
  int text_x, i, sel_lo, sel_hi;
  int orig;

  if ( NOT can_draw() )
    return;

  if ( debug )    dump( stdout, "-> DRAW_TEXT" );

  orig = set_to_glut_window();

  if ( NOT draw_text_only ) {
    if ( enabled )
      glColor3f( 1., 1., 1. );
    else
      set_to_bkgd_color();
    glDisable( GL_CULL_FACE );
    glBegin( GL_QUADS );
    glVertex2i( text_x_offset+2, 2 );     glVertex2i( w-2, 2 );
    glVertex2i( w-2, h-2 );               glVertex2i( text_x_offset+2, h-2 );
    glEnd();
  }

  /** Find where to draw the text **/

  text_x = text_x_offset + 2 + GLUI_EDITTEXT_BOXINNERMARGINX;

  /*printf( "text_x: %d      substr_width: %d     start/end: %d/%d\n",
    text_x,     substring_width( substring_start, substring_end ),
    substring_start, substring_end );
    */
  /** Find lower and upper selection bounds **/
  sel_lo = MIN(sel_start, sel_end );
  sel_hi = MAX(sel_start, sel_end );

  int sel_x_start, sel_x_end, delta;

  /** Draw selection area dark **/
  if ( sel_start != sel_end ) {
    sel_x_start = text_x;
    sel_x_end   = text_x;
    for( i=substring_start; i<=substring_end; i++ ) {
      delta = char_width( text[i] );

      if ( i < sel_lo ) {
	sel_x_start += delta;
	sel_x_end   += delta;
      }
      else if ( i < sel_hi ) {
	sel_x_end   += delta;
      }
    }

    glColor3f( 0.0f, 0.0f, .6f );
    glBegin( GL_QUADS );
    glVertex2i( sel_x_start, 2 );    glVertex2i( sel_x_end, 2 );
    glVertex2i( sel_x_end, h-2 );    glVertex2i( sel_x_start, h-2 );
    glEnd();
  }
   

  if ( sel_start == sel_end ) {   /* No current selection */
    if ( enabled )
      glColor3b( 0, 0, 0 );
    else
      glColor3b( 32, 32, 32 );
      
    glRasterPos2i( text_x, 13);
    for( i=substring_start; i<=substring_end; i++ ) {
      glutBitmapCharacter( get_font(), this->text[i] );
    }
  }
  else {                          /* There is a selection */
    int x = text_x;
    for( i=substring_start; i<=substring_end; i++ ) {
      if ( IN_BOUNDS( i, sel_lo, sel_hi-1)) { /* This character is selected */
	glColor3f( 1., 1., 1. );
	glRasterPos2i( x, 13);
	glutBitmapCharacter( get_font(), this->text[i] );
      }
      else {
	glColor3f( 0., 0., 0. );
	glRasterPos2i( x, 13);
	glutBitmapCharacter( get_font(), this->text[i] );
      }
      
      x += char_width( text[i] );
    }
  }

  restore_window( orig );

  if ( debug )    dump( stdout, "<- DRAW_TEXT" );  
}


/******************************** GLUI_EditText::find_insertion_pt() *********/
/* This function returns the character numer *before which* the insertion    */
/* point goes                                                                */

int  GLUI_EditText::find_insertion_pt( int x, int y )
{
  int curr_x, i;

  /*** See if we clicked outside box ***/
  if ( x < this->x_abs + text_x_offset )
    return -1;

  /* We move from right to left, looking to see if the mouse was clicked
     to the right of the ith character */

  curr_x = this->x_abs + text_x_offset 
    + substring_width( substring_start, substring_end )
    + 2                             /* The edittext box has a 2-pixel margin */
    + GLUI_EDITTEXT_BOXINNERMARGINX;   /** plus this many pixels blank space
					 between the text and the box       **/

  /*** See if we clicked in an empty box ***/
  if ( (int)strlen( text ) == 0 ) 
    return 0;

  /** find mouse click in text **/
  for( i=substring_end; i>=substring_start; i-- ) {
    curr_x -= char_width( text[i] );

    if ( x > curr_x ) {
      /*      printf( "-> %d\n", i );              */
      
      return i+1;
    }
  }

  return 0;

  /* Well, the mouse wasn't after any of the characters...see if it's
     before the beginning of the substring */
  if ( 0 ) {
    if ( x > (x_abs + text_x_offset + 2 ) )
      return substring_start;
    
    return -1; /* Nothing found */
  }
}


/******************************** GLUI_EditText::draw_insertion_pt() *********/

void     GLUI_EditText::draw_insertion_pt( void )
{
  int curr_x, i;

  if ( NOT can_draw() )
    return;

  /*** Don't draw insertion pt if control is disabled ***/
  if ( NOT enabled )
    return;

  if ( debug )    dump( stdout, "-> DRAW_INS_PT" );

  if ( sel_start != sel_end OR insertion_pt < 0 ) {
    return;  /* Don't draw insertion point if there is a current selection */
  }

  /*    printf( "insertion pt: %d\n", insertion_pt );              */

  curr_x = this->x_abs + text_x_offset 
    + substring_width( substring_start, substring_end )
    + 2                             /* The edittext box has a 2-pixel margin */
    + GLUI_EDITTEXT_BOXINNERMARGINX;   /** plus this many pixels blank space
					 between the text and the box       **/

  for( i=substring_end; i>=insertion_pt; i-- ) {
    curr_x -= char_width( text[i] ); 
  }  

  glColor3f( 0.0, 0.0, 0.0 );
  glBegin( GL_LINE_LOOP );
  /***
    glVertex2i( curr_x, y_abs + 4 );
    glVertex2i( curr_x, y_abs + 4 );
    glVertex2i( curr_x, y_abs + h - 3 );
    glVertex2i( curr_x, y_abs + h - 3 );
    ***/
  curr_x -= x_abs;
  glVertex2i( curr_x, 0 + 4 );
  glVertex2i( curr_x, 0 + 4 );
  glVertex2i( curr_x, 0 + h - 3 );
  glVertex2i( curr_x, 0 + h - 3 );
  glEnd();

  if ( debug )    dump( stdout, "-> DRAW_INS_PT" );
}



/******************************** GLUI_EditText::substring_width() *********/

int  GLUI_EditText::substring_width( int start, int end )
{
  int i, width;

  width = 0;

  for( i=start; i<=end; i++ )
    width += char_width( text[i] ); 

  return width;
}
 

/***************************** GLUI_EditText::update_and_draw_text() ********/

void   GLUI_EditText::update_and_draw_text( void )
{
  if ( NOT can_draw() )
    return;

  update_substring_bounds();
  /*  printf( "ss: %d/%d\n", substring_start, substring_end );                  */

  translate_and_draw_front();
}


/********************************* GLUI_EditText::special_handler() **********/

int    GLUI_EditText::special_handler( int key,int modifiers )
{
  if ( NOT glui )
    return false;
  
  if ( debug )
    printf( "SPECIAL:%d - mod:%d   subs:%d/%d  ins:%d  sel:%d/%d\n", 
	    key, modifiers, substring_start, substring_end,insertion_pt,
	    sel_start, sel_end );	 

  if ( key == GLUT_KEY_LEFT ) {
    if ( (modifiers & GLUT_ACTIVE_CTRL) != 0 ) {
      insertion_pt = find_word_break( insertion_pt, -1 );
    }
    else {
      insertion_pt--;
    }
  }
  else if ( key == GLUT_KEY_RIGHT ) {
    if ( (modifiers & GLUT_ACTIVE_CTRL) != 0 ) {
      insertion_pt = find_word_break( insertion_pt, +1 );
    }
    else {
      insertion_pt++;
    }
  }
  else if ( key == GLUT_KEY_HOME ) {
    insertion_pt = 0;
  }
  else if ( key == GLUT_KEY_END ) {
    insertion_pt = (int)strlen( text );
  }

  /*** Update selection if shift key is down ***/
  if ( (modifiers & GLUT_ACTIVE_SHIFT ) != 0 )
    sel_end = insertion_pt;
  else 
    sel_start = sel_end = insertion_pt;
  

  CLAMP( insertion_pt, 0, (int)strlen( text )); /* Make sure insertion_pt 
						   is in bounds */
  CLAMP( sel_start, 0, (int)strlen( text )); /* Make sure insertion_pt 
						is in bounds */
  CLAMP( sel_end, 0, (int)strlen( text )); /* Make sure insertion_pt 
					      is in bounds */
					      
  /******** Now redraw text ***********/
  if ( can_draw())
    update_and_draw_text();

  return true;
}


/****************************** GLUI_EditText::find_word_break() **********/
/* It looks either left or right (depending on value of 'direction'       */
/* for the beginning of the next 'word', where word are characters        */
/* separated by one of the following tokens:  " :-.,"                     */
/* If there is no next word in the specified direction, this returns      */
/* the beginning of 'text', or the very end.                              */

int    GLUI_EditText::find_word_break( int start, int direction )
{
  int    i, j;
  char   *breaks = " :-.,";
  int     num_break_chars = (int)strlen(breaks), text_len = (int)strlen(text);
  int     new_pt;

  /** If we're moving left, we have to start two back, in case we're either
    already at the beginning of a word, or on a separating token.  
    Otherwise, this function would just return the word we're already at **/
  if ( direction == -1 ) {
    start -= 2;
  }

  /***** Iterate over text in the specified direction *****/
  for ( i=start; i >= 0 AND i < text_len; i += direction ) {

    /** For each character in text, iterate over list of separating tokens **/
    for( j=0; j<num_break_chars; j++ ) {
      if ( text[i] == breaks[j] ) {

	/** character 'i' is a separating token, so we return i+1 **/
	new_pt = i + 1;

	CLAMP( new_pt, 0, text_len );

	return new_pt;
      }
    }
  }

  if ( direction > 0 )  /* Return the end of string */
    return text_len;
  else                  /* Return the beginning of the text */
    return 0;
}


/********************************** GLUI_EditText::clear_substring() ********/

void   GLUI_EditText::clear_substring( int start, int end )
{
  int i, leftover;

  /*
    printf( "clearing: %d-%d   '", start,end);
    for(i=start;i<end;i++ )
    putchar(text[i]);
    printf( "'\n" ); flushout;
    */
  /*** See if we're deleting a period in a float data-type box ***/
  if ( data_type == GLUI_EDITTEXT_FLOAT ) {
    for( i=start; i<end; i++ )
      if ( text[i] == '.' )
	num_periods = 0;
  }
  
  /*** Shift over string ***/
  leftover = (int)strlen(text) - (end);
  
  /*  printf( "leftover: %d     - ", leftover );              */

  for( i=0; i<leftover+1; i++ )
    text[start+i] = text[end+i];    

  /*  printf( "final string: '%s'\n", text );              */
}



/************************************ GLUI_EditText::update_size() **********/

void   GLUI_EditText::update_size( void )
{
  int text_size, delta;

  if ( NOT glui )
    return;

  text_size = string_width( name );

  delta = 0;
  if ( text_x_offset < text_size +2 )
    delta = text_size+2-text_x_offset;

  text_x_offset += delta;
  /*  w += delta;              */

  if ( data_type == GLUI_EDITTEXT_TEXT OR 
       data_type == GLUI_EDITTEXT_FLOAT) {
    if ( w < text_x_offset + GLUI_EDITTEXT_MIN_TEXT_WIDTH )
      w = text_x_offset + GLUI_EDITTEXT_MIN_TEXT_WIDTH;
  }
  else if ( data_type == GLUI_EDITTEXT_INT ) {
    if ( w < text_x_offset + GLUI_EDITTEXT_MIN_INT_WIDTH )
      w = text_x_offset + GLUI_EDITTEXT_MIN_INT_WIDTH;
  }
}


/****************************** GLUI_EditText::set_text() **********/

void    GLUI_EditText::set_text( char *new_text )
{
  strncpy(text,new_text,sizeof(GLUI_String));
  substring_start = 0;
  substring_end   = (int)strlen( text ) - 1;
  insertion_pt    = -1;
  sel_start       = 0;
  sel_end         = 0;

  if ( can_draw() )
    update_and_draw_text();

  /** Update the spinner, if we have one **/
  if ( spinner ) {
    spinner->float_val = this->float_val;
    spinner->int_val   = this->int_val;
  }

  /*** Now update the live variable ***/
  output_live(true);
}


/******************************* GLUI_EditText::set_float_val() ************/

void   GLUI_EditText::set_float_val( float new_val )
{
  if ( has_limits == GLUI_LIMIT_CLAMP ) {
    /*** Clamp the new value to the existing limits ***/

    CLAMP( new_val, float_low, float_high );
  } 
  else if ( has_limits == GLUI_LIMIT_WRAP ) {
    /*** Clamp the value cyclically to the limits - that is, if the
      value exceeds the max, set it the the minimum, and conversely ***/

    if ( new_val < float_low )
      new_val = float_high;
    if ( new_val > float_high )
      new_val = float_low;
  }

  float_val = new_val;
  int_val   = (int) new_val;  /* Mirror the value as an int, too */
  
  set_numeric_text();
}


/********************************** GLUI_EditText::set_int_val() ************/

void   GLUI_EditText::set_int_val( int new_val )
{
  if ( has_limits == GLUI_LIMIT_CLAMP ) {
    /*** Clamp the new value to the existing limits ***/

    CLAMP( new_val, int_low, int_high );
  }
  else if ( has_limits == GLUI_LIMIT_WRAP ) {
    /*** Clamp the value cyclically to the limits - that is, if the
      value exceeds the max, set it the the minimum, and conversely ***/

    if ( new_val < int_low )
      new_val = int_high;
    if ( new_val > int_high )
      new_val = int_low;
  }

  int_val   = new_val;
  float_val = (float) new_val;   /* We mirror the value as a float, too */

  set_numeric_text();
}


/********************************* GLUI_EditText::set_float_limits() *********/

void GLUI_EditText::set_float_limits( float low, float high, int limit_type )
{
  has_limits  = limit_type;
  float_low   = low;
  float_high  = high;
  
  if ( NOT IN_BOUNDS( float_val, float_low, float_high ))
    set_float_val( float_low );

  int_low     = (int) float_low;
  int_high    = (int) float_high;
}


/*********************************** GLUI_EditText::set_int_limits() *********/

void   GLUI_EditText::set_int_limits( int low, int high, int limit_type )
{
  has_limits  = limit_type;
  int_low     = low;
  int_high    = high;

  if ( NOT IN_BOUNDS( int_val, int_low, int_high ))
    set_int_val( int_low );

  float_low   = (float) int_low;
  float_high  = (float) int_high;
}


/************************************ GLUI_EditText::set_numeric_text() ******/

void    GLUI_EditText::set_numeric_text( void )
{
  char buf_float[200], buf_int[200];
  int  i, text_len;

  if ( data_type == GLUI_EDITTEXT_FLOAT ) {
    sprintf( buf_float, "%#g", float_val );
  
    num_periods = 0;
    text_len = (int)strlen( buf_float );
    for ( i=0; i<text_len; i++ )
      if ( buf_float[i] == '.' )
	num_periods++;

    /* Now remove trailing zeros */
    if ( num_periods > 0 ) {
      text_len = (int)strlen( buf_float );
      for ( i=text_len-1; i>0; i-- ) {
	if ( buf_float[i] == '0' AND buf_float[i-1] != '.' )
	  buf_float[i] = '\0';
	else 
	  break;
      }
    }
    
    set_text( buf_float );
  }
  else {
    sprintf( buf_int, "%d", int_val );
  
    set_text( buf_int );
  }
    
}


/*************************************** GLUI_EditText::dump() **************/

void   GLUI_EditText::dump( FILE *out, char *name )
{
  fprintf( out, 
	   "%s (edittext@%p):  ins_pt:%d  subs:%d/%d  sel:%d/%d   len:%d\n",
	   name, this, 
	   insertion_pt, substring_start, substring_end, sel_start, sel_end,
	   (int)strlen( text ));
}


/**************************************** GLUI_EditText::mouse_over() ********/

int    GLUI_EditText::mouse_over( int state, int x, int y )
{
  if ( state ) {
    /*  curr_cursor = GLUT_CURSOR_TEXT;              */
    glutSetCursor( GLUT_CURSOR_TEXT );
  }
  else {
    /*    printf( "OUT\n" );              */
    glutSetCursor( GLUT_CURSOR_LEFT_ARROW );
  }

  return true;
}
