/*
 * Copyright (C) 2011 Georgia Institute of Technology, University of Utah,
 * Weill Cornell Medical College
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This is a template header file for a user modules derived from
 * DefaultGUIModel with a custom GUI.
 */

#include <default_gui_model.h>

class HindmarshRoseModelPlugin : public DefaultGUIModel
{

  Q_OBJECT

public:
  HindmarshRoseModelPlugin(void);
  virtual ~HindmarshRoseModelPlugin(void);

  void execute(void);
  void createGUI(DefaultGUIModel::variable_t *, int);
  void customizeGUI(void);

protected:
  //Variables para el modelo

  double x, y, z, a, b, c, d, r, s, xO;

  //Variable de entorno de ejecucion
  double I, dt;
  double period;
  virtual void update(DefaultGUIModel::update_flags_t);

private:
  void initParameters();
  double getYValue(int kValue, double y, double dX, double k = 0);
  double getXValue(int kValue, double x, double dX);
  double getNextRungeKuta(double Xo, double Yo, double dX, double F(double, double, double[]), double args[]);
  void hindmarshRoseStep(double xO, double yO, double zO, double xOR, double t, double dT, double I, double a, double b, double c, double d, double r, double s);

private slots:
  // these are custom functions that can also be connected to events
  // through the Qt API. they must be implemented in plugin_template.cpp

  void aBttn_event(void);
  void bBttn_event(void);
};
