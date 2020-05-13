% Implicit kinematic constraints of
% palh3m1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [10x1]
%   Generalized joint coordinates (joint angles)
% pkin [16x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi410,phi78,phi79]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-20 17:32
% Revision: bf71b3c88fbe513355f3920541c47e1db11bd916 (2020-04-17)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = palh3m1IC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [10 1]), ...
  'palh3m1IC_kinconstr_impl_mdh_sym_varpar: qJ has to be [10x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [16 1]), ...
  'palh3m1IC_kinconstr_impl_mdh_sym_varpar: pkin has to be [16x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-20 17:27:33
% EndTime: 2020-04-20 17:27:33
% DurationCPUTime: 0.05s
% Computational Cost: add. (39->23), mult. (14->14), div. (0->0), fcn. (14->14), ass. (0->5)
t4 = -qJ(7) + pkin(15);
t3 = pkin(16) + qJ(7) + qJ(2);
t2 = -qJ(8) + t4;
t1 = qJ(3) + qJ(4) + pkin(14);
t5 = [pkin(2) * cos(t3) + cos(qJ(2)) * pkin(1) + pkin(12) - cos(qJ(6)) * pkin(5) + pkin(6); pkin(2) * sin(t3) + sin(qJ(2)) * pkin(1) - sin(qJ(6)) * pkin(5) - pkin(13); -cos(t2) * pkin(7) + pkin(3) * cos(t4) + pkin(9) * cos(t1) + cos(qJ(3)) * pkin(4); sin(t2) * pkin(7) - pkin(3) * sin(t4) + pkin(9) * sin(t1) + sin(qJ(3)) * pkin(4); -qJ(9) + qJ(6) + pi - t3; qJ(10) + pi + t1 + t2;];
h = t5;
