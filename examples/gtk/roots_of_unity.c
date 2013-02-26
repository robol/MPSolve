/**
 * Simple GTK program to demonstrate possible strategies
 * for MPSolve integration inside a graphical application.
 *
 * It uses MPSolve to solve the polynomial x^n - 1 with the n
 * specfied by the user and then displays the roots on a plot.
 *
 * Author: Leonardo Robol <robol@mail.dm.unipi.it>
 *
 * Copyright (C) 2001-2013, Dipartimento di Matematica "L. Tonelli", Pisa.
 * License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
 */

#include <mps/mps.h>
#include <gtk/gtk.h>
#include <cairo/cairo.h>

#if GTK_MAJOR_VERSION < 3

#ifndef gtk_widget_get_allocated_width
#define gtk_widget_get_allocated_width(widget) (widget->allocation.width)
#endif

#ifndef gtk_widget_get_allocated_height
#define gtk_widget_get_allocated_height(widget) (widget->allocation.height)
#endif

#endif

cplx_t * points = NULL;

int degree = 0;

GtkWidget * drawing_area = NULL;

void on_drawing_area_draw (GtkWidget* , cairo_t*);

#if GTK_MAJOR_VERSION < 3
void
on_expose_event (GtkWidget * widget, GdkEvent * event)
{
  cairo_t * cr = gdk_cairo_create (widget->window);
  on_drawing_area_draw (widget, cr);
}
#endif

void
on_drawing_area_draw (GtkWidget * widget,
                      cairo_t * cr)
{
  int width, height;

  width = gtk_widget_get_allocated_width (widget);
  height = gtk_widget_get_allocated_height (widget);

  /* Draw background */
  cairo_rectangle (cr, 0, 0, width, height);
  cairo_set_source_rgb (cr, 1, 1, 1);
  cairo_fill (cr);

  /* Draw x and y axes */
  cairo_set_line_width (cr, 1);
  cairo_set_source_rgb (cr, 0, 0, 0);

  cairo_move_to (cr, 6, 0.5 * height);
  cairo_line_to (cr, width - 6, 0.5 * height);

  cairo_move_to (cr, 0.5 * width, height - 6);
  cairo_line_to (cr, 0.5 * width, 6);

  cairo_stroke (cr);

  /* Draw points if present */
  if (points)
    {
      int i;
      double x, y;
      cairo_set_source_rgb (cr, 0.9, 0.1, 0.1);
#define PADDING (24)

      /* Check if we can draw in here. */
      if (width < 2 * PADDING || height < 2 * PADDING)
        return;

      for (i = 0; i < degree; i++)
        {
          x = cplx_Re (points[i]) * (0.5 * width - PADDING) + width / 2;
          y = -cplx_Im (points[i]) * (0.5 * height - PADDING) + height / 2;
          cairo_arc (cr, x, y, 1.3, 0, 6.29);
          cairo_fill (cr);
        }
#undef PADDING
    }
}

gboolean 
update_drawing_area (void * user_data)
{
  gtk_widget_queue_draw (drawing_area);
  return false;
}

void
on_polynomial_solved (mps_context * s, GtkButton * button)
{
  points = NULL;
  mps_context_get_roots_d (s, &points, NULL);
  degree = mps_context_get_degree (s);

  /* Call the update function from the right thread! */
  g_idle_add (update_drawing_area, NULL);

  gtk_widget_set_sensitive (GTK_WIDGET (button), true);
  mps_free_data (s);
}

void
on_solve_button_clicked (GtkButton * button, GtkSpinButton * spin_button)
{
  int degree = gtk_spin_button_get_value_as_int (spin_button);
  mps_context * s = mps_context_new ();
  mps_monomial_poly * p = mps_monomial_poly_new (s, degree);

  /* Please note that here a mutex will be required to avoid
   * freeing data while it's used to render the view. But it
   * is not added in this example to avoid complications. */
  if (points)
    {
      cplx_vfree (points);
      points = NULL;
    }

  gtk_widget_set_sensitive (GTK_WIDGET (button), false);

  /* Setting the coefficients of the input polynomial */
  mps_monomial_poly_set_coefficient_d (s, p, 0, -1, 0);
  mps_monomial_poly_set_coefficient_d (s, p, degree, 1, 0);

  /* Asking MPSolve to solve it asynchronously */
  mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));
  mps_mpsolve_async (s, (mps_callback) on_polynomial_solved, button);
}

int
main (int argc, char *argv[])
{
  /* Init the GTK environment */
  gtk_init (&argc, &argv);

  /* The window we are using for this simple example, where all the widgets will
   * be packed. */
  GtkWidget * window = gtk_window_new (GTK_WINDOW_TOPLEVEL);

  /* We will use this drawing area to draw the computed roots. */
  drawing_area = gtk_drawing_area_new ();

  /* This spin button will be used to change the degree of the computed polynomials. */
  GtkWidget * spin_button = gtk_spin_button_new_with_range (1.0, 10000.0, 1.0);
  GtkWidget * degree_label = gtk_label_new ("Degree:");
  GtkWidget * solve_button = gtk_button_new_with_label ("Solve");

#if GTK_MAJOR_VERSION < 3
  GtkWidget * degree_box = gtk_hbox_new (FALSE, 0);
#else
  GtkWidget * degree_box = gtk_box_new (GTK_ORIENTATION_HORIZONTAL, 0);
#endif

  gtk_box_pack_start (GTK_BOX (degree_box), GTK_WIDGET (degree_label), false, true, 6);
  gtk_box_pack_start (GTK_BOX (degree_box), GTK_WIDGET (spin_button), true, true, 6);
  gtk_box_pack_start (GTK_BOX (degree_box), solve_button, false, true, 6);

#if GTK_MAJOR_VERSION < 3
  GtkWidget * main_box = gtk_vbox_new (FALSE, 6);
#else
  GtkWidget * main_box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 6);
#endif
  gtk_box_pack_start (GTK_BOX (main_box), GTK_WIDGET (degree_box), false, true, 6);
  gtk_box_pack_start (GTK_BOX (main_box), GTK_WIDGET (drawing_area), true, true, 0);

  gtk_container_add (GTK_CONTAINER (window), GTK_WIDGET (main_box));
  gtk_widget_show_all (GTK_WIDGET (window));

  gtk_window_set_title (GTK_WINDOW (window), "MPSolve's roots-of-unity finder");
  gtk_widget_set_size_request (window, 420, 280);

  /* Connecting standard callbacks and the one to handle polynomial solving */
  g_signal_connect (window, "destroy", gtk_main_quit, NULL);
  g_signal_connect (solve_button, "clicked", G_CALLBACK (on_solve_button_clicked), spin_button);

#if GTK_MAJOR_VERSION < 3
  g_signal_connect (drawing_area, "expose-event", G_CALLBACK (on_expose_event), NULL);
#else
  g_signal_connect (drawing_area, "draw", G_CALLBACK (on_drawing_area_draw), NULL);
#endif

  gtk_main ();

  return EXIT_SUCCESS;
}
