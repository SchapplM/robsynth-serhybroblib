% Calculate Gravitation load on the joints for
% fourbar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:15
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbar1IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1IC_gravloadJ_floatb_twist_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:15:50
% EndTime: 2020-04-24 20:15:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (29->16), mult. (39->28), div. (4->3), fcn. (13->9), ass. (0->4)
t11 = pkin(2) * m(3);
t8 = -qJ(3) + qJ(1);
t7 = qJ(1) + qJ(2);
t1 = [(-g(2) * t11 + m(2) * (-rSges(2,1) * g(2) + rSges(2,2) * g(1))) * cos(qJ(1)) + sin(qJ(1)) * (g(1) * t11 + m(2) * (rSges(2,1) * g(1) + rSges(2,2) * g(2))) + (sin(qJ(2)) / pkin(4) * ((-rSges(4,1) * g(2) + rSges(4,2) * g(1)) * cos(qJ(3)) + (rSges(4,1) * g(1) + rSges(4,2) * g(2)) * sin(qJ(3))) * m(4) - sin(t8) / pkin(3) * ((-rSges(3,1) * g(2) + rSges(3,2) * g(1)) * cos(t7) + (rSges(3,1) * g(1) + rSges(3,2) * g(2)) * sin(t7)) * m(3)) / sin(qJ(2) + t8) * pkin(2);];
taug = t1(:);
