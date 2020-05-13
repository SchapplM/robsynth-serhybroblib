% Calculate Gravitation load on the joints for
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
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
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fivebar1IC_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1IC_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:10
% EndTime: 2020-04-27 06:19:10
% DurationCPUTime: 0.12s
% Computational Cost: add. (74->37), mult. (104->56), div. (8->4), fcn. (32->17), ass. (0->14)
t16 = qJ(1) + qJ(2);
t32 = ((-rSges(3,1) * g(2) + rSges(3,2) * g(1)) * cos(t16) + (rSges(3,1) * g(1) + rSges(3,2) * g(2)) * sin(t16)) * m(3);
t15 = qJ(3) + qJ(4);
t30 = ((-rSges(5,1) * g(2) + rSges(5,2) * g(1)) * cos(t15) + (rSges(5,1) * g(1) + rSges(5,2) * g(2)) * sin(t15)) * m(5);
t29 = pkin(2) * m(3);
t28 = pkin(3) * m(5);
t27 = cos(0.2e1 * t15) - cos(0.2e1 * t16);
t22 = 0.1e1 / pkin(5) * t32;
t21 = 0.1e1 / pkin(4) * t30;
t18 = 2 * qJ(1);
t17 = 2 * qJ(3);
t4 = -0.1e1 / sin(t15 - t16);
t3 = 0.1e1 / t27;
t1 = [(-g(2) * t29 + m(2) * (-rSges(2,1) * g(2) + rSges(2,2) * g(1))) * cos(qJ(1)) + sin(qJ(1)) * (g(1) * t29 + m(2) * (rSges(2,1) * g(1) + rSges(2,2) * g(2))) + (t27 * pkin(5) + (-cos((qJ(2) + 2 * qJ(4) + t17)) + cos((qJ(2) + t18))) * pkin(2)) * t3 * t22 + pkin(2) * sin(qJ(2)) * t4 * t21 - t32; -pkin(3) * sin(qJ(4)) * t4 * t22 + (-g(2) * t28 + m(4) * (-rSges(4,1) * g(2) + rSges(4,2) * g(1))) * cos(qJ(3)) + sin(qJ(3)) * (g(1) * t28 + m(4) * (rSges(4,1) * g(1) + rSges(4,2) * g(2))) - (t27 * pkin(4) + (cos((qJ(4) + t17)) - cos((qJ(4) + t18 + 2 * qJ(2)))) * pkin(3)) * t3 * t21 + t30;];
taug = t1(:);
