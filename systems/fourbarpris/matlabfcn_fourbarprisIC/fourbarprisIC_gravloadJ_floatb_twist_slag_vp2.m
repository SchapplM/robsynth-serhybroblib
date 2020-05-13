% Calculate Gravitation load on the joints for
% fourbarprisIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% m [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [1x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = fourbarprisIC_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisIC_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:33
% EndTime: 2020-05-07 09:59:33
% DurationCPUTime: 0.07s
% Computational Cost: add. (19->16), mult. (33->23), div. (4->4), fcn. (26->4), ass. (0->7)
t6 = mrSges(3,1) - mrSges(2,2);
t5 = m(3) * qJ(2) + mrSges(2,1) + mrSges(3,3);
t4 = cos(qJ(1));
t3 = cos(qJ(3));
t2 = sin(qJ(1));
t1 = sin(qJ(3));
t7 = [(-t3 * t1 - t4 * t2) / (pkin(3) + qJ(2)) / (-t3 ^ 2 + t4 ^ 2) * ((t6 * t2 + t5 * t4) * g(2) + (-t5 * t2 + t6 * t4) * g(1)) - 0.1e1 / pkin(2) / (t1 * t4 - t2 * t3) * (-g(1) * (t1 * mrSges(4,1) + t3 * mrSges(4,2)) - g(2) * (-t3 * mrSges(4,1) + t1 * mrSges(4,2))) + (g(1) * t4 + g(2) * t2) * m(3);];
taug = t7(:);
