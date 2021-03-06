% Calculate potential energy for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = fivebar1OL_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'fivebar1OL_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1OL_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energypot_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_energypot_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1OL_energypot_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:12:55
% EndTime: 2020-04-27 06:12:55
% DurationCPUTime: 0.11s
% Computational Cost: add. (47->40), mult. (56->46), div. (0->0), fcn. (12->8), ass. (0->9)
t12 = g(1) * r_base(1) + g(2) * r_base(2);
t11 = pkin(1) * g(1) + g(3) * r_base(3) + t12;
t9 = cos(qJ(1));
t8 = cos(qJ(3));
t7 = sin(qJ(1));
t6 = sin(qJ(3));
t5 = qJ(1) + qJ(2);
t4 = qJ(3) + qJ(4);
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * ((rSges(2,1) * g(1) + rSges(2,2) * g(2)) * t9 + (rSges(2,1) * g(2) - rSges(2,2) * g(1)) * t7 + g(3) * (r_base(3) + rSges(2,3)) + t12) - ((-rSges(3,1) * g(1) - rSges(3,2) * g(2)) * cos(t5) + (-rSges(3,1) * g(2) + rSges(3,2) * g(1)) * sin(t5) + g(3) * (r_base(3) + rSges(3,3)) + (t9 * g(1) + t7 * g(2)) * pkin(2) + t12) * m(3) - ((rSges(4,1) * g(1) + rSges(4,2) * g(2)) * t8 + (rSges(4,1) * g(2) - rSges(4,2) * g(1)) * t6 + g(3) * rSges(4,3) + t11) * m(4) - m(5) * ((rSges(5,1) * g(1) + rSges(5,2) * g(2)) * cos(t4) + (rSges(5,1) * g(2) - rSges(5,2) * g(1)) * sin(t4) + g(3) * rSges(5,3) + (t8 * g(1) + t6 * g(2)) * pkin(3) + t11);
U = t1;
