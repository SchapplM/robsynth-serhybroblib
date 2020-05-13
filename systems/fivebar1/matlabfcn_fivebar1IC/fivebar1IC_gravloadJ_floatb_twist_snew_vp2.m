% Calculate Gravitation load with newton euler on the joints for
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
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [2x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fivebar1IC_gravloadJ_floatb_twist_snew_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1IC_gravloadJ_floatb_twist_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1IC_gravloadJ_floatb_twist_snew_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_snew_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:11
% EndTime: 2020-04-27 06:19:11
% DurationCPUTime: 0.09s
% Computational Cost: add. (86->45), mult. (124->62), div. (8->4), fcn. (88->15), ass. (0->26)
t61 = sin(qJ(3));
t65 = cos(qJ(3));
t53 = t61 * g(1) - t65 * g(2);
t55 = -t65 * g(1) - t61 * g(2);
t60 = sin(qJ(4));
t64 = cos(qJ(4));
t76 = (mrSges(5,1) * (t64 * t53 - t60 * t55) - mrSges(5,2) * (t60 * t53 + t64 * t55)) / pkin(4);
t63 = sin(qJ(1));
t67 = cos(qJ(1));
t54 = t63 * g(1) - t67 * g(2);
t56 = -t67 * g(1) - t63 * g(2);
t62 = sin(qJ(2));
t66 = cos(qJ(2));
t75 = (mrSges(3,1) * (-t66 * t54 + t62 * t56) - mrSges(3,2) * (-t62 * t54 - t66 * t56)) / pkin(5);
t72 = qJ(4) + qJ(3);
t73 = qJ(1) + qJ(2);
t74 = cos(0.2e1 * t72) - cos(0.2e1 * t73);
t69 = 0.2e1 * qJ(1);
t68 = 0.2e1 * qJ(3);
t57 = -0.1e1 / sin(t72 - t73);
t52 = t62 * mrSges(3,1) + t66 * mrSges(3,2) - mrSges(2,2);
t51 = t60 * mrSges(5,1) + t64 * mrSges(5,2) + mrSges(4,2);
t50 = pkin(2) * m(3) - t66 * mrSges(3,1) + t62 * mrSges(3,2) + mrSges(2,1);
t49 = pkin(3) * m(5) + t64 * mrSges(5,1) - t60 * mrSges(5,2) + mrSges(4,1);
t48 = 0.1e1 / t74;
t1 = [-(-t63 * t50 + t52 * t67) * g(1) - (t50 * t67 + t63 * t52) * g(2) - (t74 * pkin(5) + (-cos(qJ(2) + 0.2e1 * qJ(4) + t68) + cos(qJ(2) + t69)) * pkin(2)) * t48 * t75 + pkin(2) * t62 * t57 * t76; pkin(3) * t60 * t57 * t75 - (-t61 * t49 - t51 * t65) * g(1) - (t49 * t65 - t51 * t61) * g(2) - (t74 * pkin(4) + (cos(qJ(4) + t68) - cos(qJ(4) + t69 + 0.2e1 * qJ(2))) * pkin(3)) * t48 * t76;];
taug = t1(:);
